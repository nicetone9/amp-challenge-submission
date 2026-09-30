"""Cold-clone the public repository and run default frozen-pool filtering twice."""
import argparse
import json
from pathlib import Path
import runpy
import subprocess

from amp_submission.io import digest, read_fasta, verify_assets
from amp_submission.prepare_reference import save_json
from reproduce_from_model import REPOSITORY, cold_environment, runtime_command


def reproduce(root, revision):
    if not revision or len(revision) != 40 or set(revision) - set("0123456789abcdef"):
        raise ValueError("Use a full immutable Git commit SHA")
    if root.exists():
        raise ValueError("Use a new verification directory")
    root.mkdir(parents=True)
    reports = []
    for name in ("run1", "run2"):
        target = root / name
        env = dict(cold_environment(target), UV_FROZEN="1", PYTHONHASHSEED="42",
                   OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4",
                   TOKENIZERS_PARALLELISM="false", AMP_TRACK_ONLINE="1")
        save_json(root / "status.json", {"status": "RUNNING", "run": name})
        with (root / (name + ".log")).open("w") as log:
            for command in (
                ["git", "clone", "--no-checkout", REPOSITORY, str(target)],
                ["git", "-C", str(target), "checkout", "--detach", revision],
                [*runtime_command(target), "uv", "sync", "--frozen"],
                [*runtime_command(target), "uv", "run", "generate"],
            ):
                result = subprocess.run(command, cwd=target if target.exists() else root,
                                        env=env, stdout=log, stderr=subprocess.STDOUT)
                if result.returncode:
                    save_json(root / "status.json", {"status": "FAILED", "run": name,
                                                     "command": command, "exit_code": result.returncode})
                    raise RuntimeError(name + " failed; dependent run stopped")
        actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=target, text=True).strip()
        if actual != revision:
            raise ValueError("Clone revision changed")
        report = json.loads((target / "generate/generation.json").read_text())
        if report["fresh_model_sampling"] or not report["candidate_cache_read"]:
            raise ValueError("Default command unexpectedly resampled models")
        if (report["pool_size"], report["n_sequences"], report["top_k"]) != (600000, 50000, 100):
            raise ValueError("Nondefault quantities cannot certify this pool")
        reports.append(report)
    for key in ("seed", "pool_size", "n_sequences", "top_k", "motif_quota", "thresholds",
                "quality_metrics", "descriptive_only", "threshold", "pool_manifest_sha256",
                "reference_sha256", "code_sha256", "locks", "raw_sha256", "environment"):
        if reports[0][key] != reports[1][key]:
            raise ValueError("Two-run mismatch: " + key)
    hashes = {}
    provenance = json.loads((root / "run1/frozen-pool/provenance.json").read_text())
    validators = runpy.run_path(str(root / "run1/scripts/verify_submission.py"))
    for filename in ("library.fasta", "top.fasta"):
        first, second = (root / name / "generate" / filename for name in ("run1", "run2"))
        if first.read_bytes() != second.read_bytes():
            raise ValueError("Byte identity failed: " + filename)
        hashes[filename] = digest(first)
        if hashes[filename] != provenance["source_fasta_sha256"][filename.removesuffix(".fasta")]:
            raise ValueError("Frozen filtering changed the source pool finalization selection")
    for name in ("run1", "run2"):
        target = root / name
        verify_assets(target / "frozen-pool")
        output = target / "generate"
        references = set(read_fasta(target / "frozen-pool/antibacterial.fasta"))
        library = validators["_verify_sequences"](output / "library.fasta")
        validators["_verify_top"](output / "top.fasta", library, 100)
        validators["_verify_no_overlap"](library, references)
        validators["_veritfy_max_simularity"](set(read_fasta(output / "top.fasta")), references)
    result = {"status": "TWO_PUBLIC_COLD_CLONES_FROZEN_FILTER_BYTE_IDENTICAL",
              "revision": revision, "command_each_run": ["uv", "run", "generate"],
              "public_inputs_downloaded_each_clone": True,
              "fresh_model_sampling": False, "frozen_source_selection_match": True,
              "pool_counts": {"vq": 300000, "dima": 300000}, "seed": 42,
              "library_count": 50000, "top_count": 100, "fasta_sha256": hashes,
              "raw_sha256": reports[0]["raw_sha256"],
              "pool_manifest_sha256": reports[0]["pool_manifest_sha256"],
              "threshold": reports[0]["threshold"], "environment": reports[0]["environment"],
              "official_template_component_checks_passed": True,
              "swanlab_run_ids": [report["swanlab_run_id"] for report in reports],
              "not_claimed": ["model-resampling identity", "organizer acceptance",
                              "six-category quality completion", "wet-lab activity"]}
    save_json(root / "verification.json", result)
    save_json(root / "status.json", {"status": "VERIFIED"})
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("work/frozen-pool-cold-reproduction"))
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    reproduce(args.root.resolve(), args.revision)
