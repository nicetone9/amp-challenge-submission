"""Finish fixed-policy scoring/selection from checksum-verified completed raw samples."""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import shutil
import time

import numpy as np
import pandas as pd
import swanlab

from amp_submission.generate import parser as generate_parser, validate
from amp_submission.io import digest, verify_assets, write_fasta
from amp_submission.model_generate import choose, load_reference, merge_pool, score_candidates
from amp_submission.prepare_reference import save_json
from freeze_candidate_pool import load_completed_raw


def finalize(source, assets, output):
    source, assets, output = map(Path, (source, assets, output))
    if output.exists():
        raise ValueError("Use a new finalization directory")
    inputs, complete, branches = load_completed_raw(source)
    verify_assets(assets)
    reference = load_reference(assets / "generation-reference")
    if (digest(assets / "manifest.json") != inputs["assets_sha256"]
            or reference["manifest_sha256"] != inputs["reference_sha256"]):
        raise ValueError("Source assets/reference changed")
    args = generate_parser().parse_args([])
    for key in ("seed", "pool_size", "threads", "n_sequences", "top_k", "motif_quota"):
        if getattr(args, key) != inputs[key]:
            raise ValueError("Recovery cannot change the original selection policy: " + key)
    cdhit = shutil.which("cd-hit")
    if cdhit is None:
        raise RuntimeError("CD-HIT is missing")
    output.mkdir(parents=True)
    shutil.copytree(source / "raw", output / "raw")
    selection_code = {str(p.relative_to(Path("src/amp_submission"))): digest(p)
                      for p in sorted(Path("src/amp_submission").rglob("*.py"))}
    config = {**inputs, "reused_completed_raw": True, "sampled_in_this_invocation": False,
              "selection_code_sha256": selection_code,
              "failed_sampling_run_id": "lduq0g39"}
    save_json(output / "inputs.json", config)
    run = swanlab.init(workspace="nicetone9", project="AMP_step2challenge", mode="online",
                       name="finalize-completed-600k", group="frozen-pool-seed42",
                       job_type="selection-recovery", config=config, log_dir=str(output / "swanlog"))
    start, step = time.monotonic(), 0
    def progress(key, value):
        nonlocal step
        event = {"stage": key, "completed": value, "step": step,
                 "wall_seconds": time.monotonic() - start}
        print(json.dumps(event), flush=True)
        save_json(output / "progress.json", event)
        swanlab.log({key: value, "cost/wall_seconds": event["wall_seconds"]}, step=step)
        step += 1
    try:
        progress("input/verified_raw_candidates", sum(len(s) for s in branches.values()))
        rows = merge_pool(branches, reference["training"])
        frame, percentiles, hard = score_candidates(rows, reference, output, progress)
        library, top, threshold, attempts, command = choose(
            frame, percentiles, hard, reference["motifs"], args, output)
        validate([r["sequence"] for r in library], [r["sequence"] for r in top],
                 reference["known"], args.n_sequences, args.top_k)
        for name, selected in (("library", library), ("top", top)):
            write_fasta(output / (name + ".fasta"), [r["sequence"] for r in selected])
            pd.DataFrame(selected).to_csv(output / (name + "-scores.csv"), index=False)
        report = {**config, "status": "RECOVERED_SAMPLED_POOL_HARD_COMPLIANT_PARTIAL_QUALITY",
                  "all_six_categories_passed": False, "threshold": threshold, "attempts": attempts,
                  "pool_unique": len(rows),
                  "source_counts": dict(Counter(r["source"] for r in library)),
                  "tier_counts": dict(Counter(r["tier"] for r in library)),
                  "motif_supported_library": sum(bool(r["short_motif_support"]) for r in library),
                  "motif_supported_top": sum(bool(r["short_motif_support"]) for r in top),
                  "fasta_sha256": {name: digest(output / (name + ".fasta")) for name in ("library", "top")},
                  "raw_sha256": complete["sha256"],
                  "environment": {"python": platform.python_version(), "numpy": np.__version__,
                                  "pandas": pd.__version__, "device": "cpu",
                                  "cd_hit_sha256": digest(cdhit)},
                  "cd_hit_command": command, "wall_seconds": time.monotonic() - start,
                  "swanlab_run_id": run.id,
                  "warning": "Raw sampling was complete; descriptors were repaired before CPU selection. "
                             "Not a successful fresh-model cold run or six-category/wet-lab validation."}
        save_json(output / "generation.json", report)
        progress("generation/library_count", len(library))
        progress("generation/top_count", len(top))
        run.finish()
        print(json.dumps({"status": report["status"], "output": str(output),
                          "threshold": threshold, "fasta_sha256": report["fasta_sha256"]}), flush=True)
    except BaseException as error:
        save_json(output / "failure.json", {"error": type(error).__name__, "detail": str(error)})
        run.finish(state="crashed", error=str(error))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--assets", type=Path, default=Path("checkpoint"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    finalize(args.source, args.assets, args.output)
