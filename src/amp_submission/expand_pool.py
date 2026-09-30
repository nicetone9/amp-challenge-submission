"""Expand the raw mixed pool in deterministic resumable batches; not qualification."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from .calibration import fingerprint, subseed
from .io import digest, read_fasta, verify_assets
from .prepare_reference import save_json

def extend_batches(output, identity, target_per_branch, generate_batch, batch_size=128, progress_callback=None):
    """Independent batch seeds make warm-cache and interrupted runs identical."""
    output = Path(output)
    binding = output / "identity.json"
    if binding.exists() and json.loads(binding.read_text()) != identity:
        raise ValueError("Pool resume identity mismatch")
    save_json(binding, identity)
    counts = {}
    # Always generate a full batch. A larger later budget cannot alter its prefix.
    for arch in ("vq", "dima"):
        count = 0
        for index in range((target_per_branch + batch_size - 1) // batch_size):
            path = output / arch / f"batch-{index:07d}.json"
            batch_identity = fingerprint({"identity": identity, "arch": arch, "index": index,
                                          "batch_size": batch_size})
            cached = path.exists()
            if cached:
                batch = json.loads(path.read_text())
                if (batch["identity"] != batch_identity or
                        batch["sequences_sha256"] != fingerprint(batch["sequences"]) or
                        len(batch["sequences"]) != batch_size):
                    raise ValueError("Corrupt pool batch: " + str(path))
            else:
                sequences = generate_batch(arch, index, batch_size,
                                           subseed(identity["seed"], "mixed-expansion", arch, index))
                if len(sequences) != batch_size:
                    raise ValueError("Sampler returned wrong batch size")
                batch = {"identity": batch_identity, "arch": arch, "index": index,
                         "sequences": sequences, "sequences_sha256": fingerprint(sequences)}
                save_json(path, batch)
            count += min(batch_size, target_per_branch - count)
            if progress_callback is not None:
                progress_callback(arch, index, count, cached)
            save_json(output / "progress.json",
                      {"identity": fingerprint(identity), "current_arch": arch,
                       "completed_full_batches": index + 1, "requested_per_branch": target_per_branch})
        counts[arch] = count
    return counts

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assets", type=Path, default=Path("checkpoint"))
    p.add_argument("--initial", type=Path, default=Path("candidates"))
    p.add_argument("--output", type=Path, default=Path("work/expanded-pool"))
    p.add_argument("--pool-size", type=int, default=300000)
    p.add_argument("--max-candidates", type=int, default=3000000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--min-length", type=int, default=8)
    p.add_argument("--max-length", type=int, default=40)
    p.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    p.add_argument("--batch-size", type=int, default=128)
    a = p.parse_args()
    if not 100000 <= a.pool_size <= a.max_candidates or a.pool_size % 2:
        raise ValueError("Even raw pool target required, between 100000 and max-candidates")
    if not 8 <= a.min_length <= a.max_length <= 40 or a.batch_size < 1:
        raise ValueError("Unsupported scoring length or batch size")
    physical_raw = 100000 + 2 * (((a.pool_size - 100000) // 2 + a.batch_size - 1) // a.batch_size) * a.batch_size
    if physical_raw > a.max_candidates:
        raise ValueError("Full deterministic batches would exceed max-candidates; reduce pool-size")
    verify_assets(a.assets)
    initial = {}
    for arch in ("vq", "dima"):
        path = a.initial / arch / "library.fasta"
        seqs = read_fasta(path)
        if len(seqs) != 50000 or len(set(seqs)) != 50000:
            raise ValueError("Expected existing unique 50000-sequence branch: " + arch)
        initial[arch] = {"path": str(path.resolve()), "sha256": digest(path), "count": len(seqs)}
    identity = {"seed": a.seed, "assets_sha256": digest(a.assets / "manifest.json"),
                "initial": initial, "batch_size": a.batch_size, "min_length": a.min_length,
                "max_length": a.max_length, "device": a.device, "pixi_lock": digest("pixi.lock"),
                "code": {name: digest(Path(__file__).parent / name) for name in
                         ("expand_pool.py", "sampling.py", "models.py")}}
    import torch
    from .sampling import OBJECTIVES, load_models, sample, seed_all
    if a.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("Allocated CUDA GPU unavailable; no silent CPU fallback")
    torch.set_num_threads(4)
    pdf = np.asarray(json.loads((a.assets / "length-pdf.json").read_text()), float)
    pdf[:a.min_length] = 0
    pdf[a.max_length+1:] = 0
    pdf /= pdf.sum()
    loaded = {}
    current = None
    def generate(arch, index, n, seed):
        nonlocal current
        if current != arch:
            loaded.clear()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            policies, reference, bundle = load_models(a.assets, arch, a.device)
            loaded.update(policies=policies, reference=reference, bundle=bundle)
            current = arch
        seed_all(seed % (2**32))
        lengths = np.random.default_rng(seed).choice(len(pdf), n, p=pdf).tolist()
        objective = OBJECTIVES[index % len(OBJECTIVES)]
        return sample(arch, loaded["policies"][objective], loaded["reference"],
                      loaded["bundle"], lengths, a.device)
    started = time.monotonic()
    import swanlab
    tracking_path = a.output / "tracking.json"
    tracking = json.loads(tracking_path.read_text()) if tracking_path.exists() else None
    run_id = fingerprint([identity, str(a.output.resolve())])[:32]
    if tracking and tracking["id"] != run_id:
        raise ValueError("SwanLab resume identity mismatch")
    run = swanlab.init(workspace="nicetone9", project="AMP_step2challenge", mode="online",
                       id=run_id, resume="must" if tracking else "allow",
                       name="mixed-pool-expansion", group="six-metrics-seed42",
                       job_type="sampling", config={**identity, "scope": "raw pool, not qualified AMPs"},
                       log_dir=str(a.output / "swanlog"))
    state = tracking or {"id": run_id, "last_step": -1}
    save_json(tracking_path, state)
    def progress(arch, index, count, cached):
        if cached or index % 16:
            return
        state["last_step"] += 1
        swanlab.log({"sampling/" + arch + "_additional_raw": count,
                     "sampling/requested_pool": a.pool_size,
                     "sampling/wall_seconds": time.monotonic() - started}, step=state["last_step"])
        save_json(tracking_path, state)
    try:
        counts = extend_batches(a.output, identity, (a.pool_size - 100000) // 2,
                                generate, a.batch_size, progress)
    except BaseException as error:
        run.finish(state="crashed", error=type(error).__name__ + ": " + str(error))
        raise
    state["last_step"] += 1
    swanlab.log({"sampling/raw_pool_complete": a.pool_size,
                 "sampling/qualified_count_known": 0}, step=state["last_step"])
    save_json(tracking_path, state)
    run.finish()
    report = {"status": "RAW_POOL_ONLY_NOT_QUALIFIED", "raw_requested": a.pool_size,
              "additional_counts": counts, "physical_raw_including_batch_tail": physical_raw,
              "max_candidates": a.max_candidates,
              "identity_sha256": fingerprint(identity), "wall_seconds": time.monotonic() - started,
              "warning": "Raw count is not unique qualified count. Six-category evaluation remains required."}
    save_json(a.output / "complete.json", report)
    print(json.dumps(report))

if __name__ == "__main__":
    main()
