"""Fresh mixed-model inference; never reads a previous candidate pool or delivery."""
import json
import os
import platform
import shutil
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import fingerprint, subseed
from .hard_select import METRICS, DESCRIPTIVE_ONLY, THRESHOLDS, quality_arrays
from .io import AA, digest, read_fasta, verify_assets, write_fasta
from .metric_scores import anchor_novelty, descriptors
from .motif_select import motif_top
from .motifs import scan
from .prepare_reference import save_json
from .selection import assign_tiers, cluster_sequences, select


def load_reference(root):
    root = Path(root)
    manifest = verify_assets(root)
    required = {"reference-scores.json", "anchors.fasta", "known-amp.fasta",
                "training-hashes.txt", "motifs.json", "provenance.json"}
    if not required <= set(manifest["sha256"]):
        raise ValueError("Incomplete frozen generation reference")
    reference = json.loads((root / "reference-scores.json").read_text())
    if not set(METRICS) <= set(reference) or set(reference) - set(METRICS) - set(DESCRIPTIVE_ONLY):
        raise ValueError("Unexpected reference metric panel")
    arrays = [np.asarray(reference[m], float) for m in METRICS]
    if len({len(a) for a in arrays}) != 1 or not len(arrays[0]) or not all(np.isfinite(a).all() for a in arrays):
        raise ValueError("Nonfinite or misaligned reference scores")
    motifs = json.loads((root / "motifs.json").read_text())["motifs"]
    return {"scores": reference, "motifs": sorted(motifs),
            "anchors": read_fasta(root / "anchors.fasta"),
            "known": read_fasta(root / "known-amp.fasta"),
            "training": set((root / "training-hashes.txt").read_text().splitlines()),
            "manifest_sha256": digest(root / "manifest.json")}


def merge_pool(branches, training):
    merged = {}
    for arch in ("vq", "dima"):
        for sequence in branches[arch]:
            if sequence in merged:
                row = merged[sequence]
                row["source"] = "|".join(sorted(set(row["source"].split("|")) | {arch}))
                row["raw_occurrences"] += 1
                continue
            valid = (8 <= len(sequence) <= 50 and set(sequence) <= set(AA))
            low_complexity = bool(sequence) and max(Counter(sequence).values()) / len(sequence) > .6
            merged[sequence] = {"sequence": sequence, "source": arch, "raw_occurrences": 1,
                                "hard_precheck": valid and not low_complexity
                                and fingerprint(sequence) not in training}
    return list(merged.values())


def sample_fresh(args, directory, progress):
    import torch
    from .sampling import OBJECTIVES, load_models, sample, seed_all
    torch.set_num_threads(args.threads)
    if not torch.cuda.is_available():
        raise RuntimeError("Default mixed generation requires an allocated CUDA GPU; no CPU fallback")
    pdf = np.asarray(json.loads((args.assets / "length-pdf.json").read_text()), float)
    pdf[:8] = 0
    pdf[51:] = 0
    if not np.isfinite(pdf).all() or np.any(pdf < 0) or pdf.sum() <= 0:
        raise ValueError("Invalid frozen length distribution")
    pdf /= pdf.sum()
    seed_all(args.seed)
    branches, hashes = {}, {}
    target, batch_size = args.pool_size // 2, 128
    for arch in ("vq", "dima"):
        policies, reference, bundle = load_models(args.assets, arch, "cuda")
        sequences = []
        for index in range((target + batch_size - 1) // batch_size):
            seed = subseed(args.seed, "mixed-expansion", arch, index)
            seed_all(seed % (2**32))
            lengths = np.random.default_rng(seed).choice(len(pdf), batch_size, p=pdf).tolist()
            batch = sample(arch, policies[OBJECTIVES[index % len(OBJECTIVES)]],
                           reference, bundle, lengths, "cuda")
            if len(batch) != batch_size:
                raise ValueError("Model returned an incomplete batch")
            save_json(directory / "raw" / arch / f"batch-{index:07d}.json",
                      {"index": index, "seed": seed, "sequences": batch,
                       "sha256": fingerprint(batch)})
            sequences.extend(batch)
            if index % 32 == 0 or len(sequences) >= target:
                progress("sampling/" + arch, min(len(sequences), target))
        branches[arch] = sequences[:target]
        hashes[arch] = fingerprint(branches[arch])
        del policies, reference, bundle
        torch.cuda.empty_cache()
    save_json(directory / "raw/complete.json",
              {"fresh_model_sampling": True, "initial_library_used": False,
               "counts": {a: len(s) for a, s in branches.items()}, "sha256": hashes,
               "physical_attempts": 2 * ((target + 127) // 128) * 128})
    return branches


def score_candidates(rows, reference, directory, progress):
    from .generate import too_similar
    sequences = [r["sequence"] for r in rows]
    if any(not 8 <= len(s) <= 50 or not set(s) <= set(AA) for s in sequences):
        raise ValueError("Model generated unsupported sequence; no truncation or silent repair")
    score_names = (*METRICS, *DESCRIPTIVE_ONLY)
    columns = {m: [] for m in score_names}
    for start in range(0, len(sequences), 2048):
        batch = sequences[start:start+2048]
        values = descriptors(batch)
        values["anchor_novelty"] = anchor_novelty(batch, reference["anchors"])
        for metric in score_names:
            columns[metric].extend(values[metric].tolist())
        progress("scoring/sequences", start + len(batch))
    frame = pd.DataFrame(rows)
    for metric in score_names:
        frame[metric] = columns[metric]
    percentiles, ranking, weakest = quality_arrays(reference["scores"], frame)
    finite = np.isfinite(percentiles).all(axis=1)
    possible = frame.hard_precheck.to_numpy(bool) & finite & (percentiles >= min(THRESHOLDS)).all(axis=1)
    # Everyone receives the same five quality scores (length and mass are descriptive only). Expensive whole-reference checks
    # are required for every selectable candidate, never inferred from percentiles.
    hard = np.zeros(len(frame), dtype=bool)
    indices = np.flatnonzero(possible)
    for j, i in enumerate(indices):
        hard[i] = not too_similar(sequences[i], reference["known"])
        if j % 2048 == 0 or j + 1 == len(indices):
            progress("hard_reference/audited", j + 1)
    frame["known_amp_ratio_pass"] = hard
    frame["known_amp_ratio_audited"] = possible
    frame["ranking"] = ranking
    frame["weakest"] = weakest
    for j, metric in enumerate(METRICS):
        frame[metric + "_percentile"] = percentiles[:, j]
    frame.to_csv(directory / "candidate-scores.csv", index=False)
    return frame, percentiles, hard


def choose(frame, percentiles, hard, motifs, args, directory):
    attempts = []
    short = sorted(m for m in motifs if 3 <= len(m) <= 4)
    for threshold in THRESHOLDS:
        indices = np.flatnonzero(hard & (percentiles >= threshold).all(axis=1))
        attempt = {"threshold": threshold, "eligible": len(indices)}
        attempts.append(attempt)
        if len(indices) < args.n_sequences:
            continue
        eligible = frame.iloc[indices].to_dict("records")
        clusters, command = cluster_sequences([r["sequence"] for r in eligible],
                                              directory / f"cluster-{threshold:.2f}")
        for row in eligible:
            row.update(cluster=clusters[row["sequence"]], passed=True,
                       motif_support=sum(m in row["sequence"] for m in motifs))
        tiers = assign_tiers(eligible, args.seed)
        capacities = Counter(r["tier"] for r in tiers)
        attempt["tier_capacity"] = dict(capacities)
        if min(capacities.get("high", 0), capacities.get("low", 0)) < args.n_sequences // 2:
            continue
        library = select(tiers, args.n_sequences)
        for row in library:
            hits = scan(row["sequence"], short)
            row["short_motif_support"] = len({h["motif"] for h in hits})
            row["short_motifs"] = ";".join(sorted({h["motif"] for h in hits}))
            row["motif_positions_json"] = json.dumps(hits, sort_keys=True, separators=(",", ":"))
        # Missing motif quota is a failure, not permission to alter the ranking policy.
        top = motif_top(library, args.top_k, args.motif_quota)
        return library, top, threshold, attempts, command
    raise ValueError("Insufficient qualified candidates within frozen threshold grid; no filler")


def run(args):
    from .generate import validate
    if args.device not in ("auto", "cuda"):
        raise ValueError("Mixed default is validated on CUDA only")
    if args.seed < 0 or args.threads < 1 or args.pool_size < 2 or args.pool_size % 2:
        raise ValueError("Require nonnegative seed, positive threads, even pool-size")
    if not 2 <= args.top_k <= args.n_sequences <= args.pool_size or args.n_sequences % 2 or args.top_k % 2:
        raise ValueError("Require even library/top sizes with 2 <= top <= library <= pool")
    if not 0 <= args.motif_quota <= args.top_k or args.motif_quota % 2:
        raise ValueError("Motif quota must be even and between zero and top-k")
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    verify_assets(args.assets)
    reference = load_reference(args.selection_reference or args.assets / "generation-reference")
    official = read_fasta(args.assets / "antibacterial.fasta")
    if not set(official) <= set(reference["known"]):
        raise ValueError("Frozen novelty reference omits official sequences")
    # A new retained staging directory every invocation: no candidate-cache reads.
    args.output.parent.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix=".amp-fresh-", dir=args.output.parent))
    config = {"seed": args.seed, "pool_size": args.pool_size, "batch_size": 128,
              "threads": args.threads, "n_sequences": args.n_sequences, "top_k": args.top_k,
              "motif_quota": args.motif_quota, "thresholds": list(THRESHOLDS),
              "quality_metrics": list(METRICS), "descriptive_only": list(DESCRIPTIVE_ONLY),
              "length_hard_bounds": [8, 50],
              "assets_sha256": digest(args.assets / "manifest.json"),
              "reference_sha256": reference["manifest_sha256"],
              "fresh_model_sampling": True, "candidate_cache_read": False,
              "code_sha256": {str(p.relative_to(Path(__file__).parent)): digest(p)
                              for p in sorted(Path(__file__).parent.rglob("*.py"))},
              "locks": {name: digest(name) for name in ("pixi.lock", "uv.lock")}}
    save_json(directory / "inputs.json", config)
    run_handle = None
    if args.track:
        import swanlab
        run_handle = swanlab.init(workspace="nicetone9", project="AMP_step2challenge",
                                 mode="online", name="fresh-model-generate",
                                 group="from-model-twice-seed42", job_type="reproducibility",
                                 config=config, log_dir=str(directory / "swanlog"))
    start, step = time.monotonic(), 0
    def progress(key, value):
        nonlocal step
        event = {"stage": key, "completed": value, "step": step,
                 "wall_seconds": time.monotonic() - start}
        print(json.dumps(event), flush=True)
        save_json(directory / "progress.json", event)
        if run_handle is not None:
            swanlab.log({key: value, "cost/wall_seconds": event["wall_seconds"]}, step=step)
        step += 1
    try:
        import torch
        branches = sample_fresh(args, directory, progress)
        rows = merge_pool(branches, reference["training"])
        frame, percentiles, hard = score_candidates(rows, reference, directory, progress)
        library, top, threshold, attempts, command = choose(
            frame, percentiles, hard, reference["motifs"], args, directory)
        seqs, top_seqs = [r["sequence"] for r in library], [r["sequence"] for r in top]
        validate(seqs, top_seqs, reference["known"], args.n_sequences, args.top_k)
        for name, selected in (("library", library), ("top", top)):
            write_fasta(directory / (name + ".fasta"), [r["sequence"] for r in selected])
            pd.DataFrame(selected).to_csv(directory / (name + "-scores.csv"), index=False)
        report = {**config, "status": "FRESH_MODEL_HARD_COMPLIANT_PARTIAL_QUALITY",
                  "all_six_categories_passed": False, "threshold": threshold,
                  "attempts": attempts, "pool_unique": len(rows),
                  "source_counts": dict(Counter(r["source"] for r in library)),
                  "tier_counts": dict(Counter(r["tier"] for r in library)),
                  "motif_supported_library": sum(bool(r["short_motif_support"]) for r in library),
                  "motif_supported_top": sum(bool(r["short_motif_support"]) for r in top),
                  "fasta_sha256": {name: digest(directory / (name + ".fasta")) for name in ("library", "top")},
                  "raw_sha256": {a: fingerprint(s) for a, s in branches.items()},
                  "environment": {"python": platform.python_version(), "torch": torch.__version__,
                                  "cuda": torch.version.cuda, "gpu": torch.cuda.get_device_name(0)},
                  "cd_hit_command": command, "wall_seconds": time.monotonic() - start,
                  "artifact_directory": str(directory.resolve()),
                  "swanlab_run_id": run_handle.id if run_handle else None,
                  "warning": "Five-metric internal ranking; length and mass are descriptive only; not official score or measured activity. "
                             "A single run is not a two-run reproducibility certificate."}
        save_json(directory / "generation.json", report)
        progress("generation/library_count", len(library))
        progress("generation/top_count", len(top))
        if run_handle:
            run_handle.finish()
        args.output.mkdir(parents=True, exist_ok=True)
        for name in ("library.fasta", "top.fasta", "library-scores.csv", "top-scores.csv", "generation.json"):
            temporary = args.output / (name + ".partial")
            shutil.copyfile(directory / name, temporary)
            os.replace(temporary, args.output / name)
        print(json.dumps(report), flush=True)
        return report
    except BaseException as error:
        save_json(directory / "failure.json", {"error": type(error).__name__, "detail": str(error)})
        if run_handle:
            run_handle.finish(state="crashed", error=str(error))
        raise
