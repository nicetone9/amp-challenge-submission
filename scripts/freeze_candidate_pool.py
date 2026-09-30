"""Freeze one completed VQ/DiMA sampling run for public deterministic filtering."""
import argparse
import gzip
import json
import shutil
import tarfile
from pathlib import Path

from amp_submission.calibration import fingerprint, subseed
from amp_submission.frozen_pool import load_pool
from amp_submission.io import digest, verify_assets, write_fasta
from amp_submission.model_generate import load_reference
from amp_submission.prepare_reference import save_json


def load_completed_raw(source):
    source = Path(source)
    complete = json.loads((source / "raw/complete.json").read_text())
    inputs = json.loads((source / "inputs.json").read_text())
    if not complete["fresh_model_sampling"] or complete["counts"] != {"vq": 300000, "dima": 300000}:
        raise ValueError("Require the completed 600,000-attempt mixed pool")
    branches = {}
    for arch in ("vq", "dima"):
        sequences = []
        paths = sorted((source / "raw" / arch).glob("batch-*.json"))
        if len(paths) != (complete["counts"][arch] + 127) // 128:
            raise ValueError("Incomplete raw model batches")
        for index, path in enumerate(paths):
            batch = json.loads(path.read_text())
            if (batch["index"] != index or len(batch["sequences"]) != 128
                    or batch["seed"] != subseed(inputs["seed"], "mixed-expansion", arch, index)
                    or batch["sha256"] != fingerprint(batch["sequences"])):
                raise ValueError("Changed raw model batch: " + str(path))
            sequences.extend(batch["sequences"])
        branches[arch] = sequences[:complete["counts"][arch]]
        if fingerprint(branches[arch]) != complete["sha256"][arch]:
            raise ValueError("Raw branch checksum mismatch")
    return inputs, complete, branches


def freeze(source, assets, output, archive):
    source, assets, output, archive = map(Path, (source, assets, output, archive))
    if output.exists() or archive.exists():
        raise ValueError("Use new frozen-pool and archive paths")
    inputs, complete, branches = load_completed_raw(source)
    report = json.loads((source / "generation.json").read_text())
    reference_root = assets / "generation-reference"
    reference_manifest = verify_assets(reference_root)
    if digest(reference_root / "manifest.json") != inputs["reference_sha256"]:
        raise ValueError("Source generation reference changed")
    verify_assets(assets)
    if digest(assets / "manifest.json") != inputs["assets_sha256"]:
        raise ValueError("Source model assets changed")
    output.mkdir(parents=True)
    for arch, sequences in branches.items():
        write_fasta(output / (arch + ".fasta"), sequences)
    shutil.copyfile(source / "candidate-scores.csv", output / "candidate-scores.csv")
    shutil.copyfile(assets / "antibacterial.fasta", output / "antibacterial.fasta")
    for filename in ("manifest.json", *reference_manifest["sha256"]):
        target = output / "reference" / filename
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(reference_root / filename, target)
    provenance = {"sampling_seed": inputs["seed"], "counts": complete["counts"],
                  "raw_sha256": complete["sha256"], "physical_attempts": complete["physical_attempts"],
                  "assets_sha256": inputs["assets_sha256"], "reference_sha256": inputs["reference_sha256"],
                  "sampling_code_sha256": inputs["code_sha256"],
                  "selection_code_sha256": report.get("selection_code_sha256", inputs["code_sha256"]),
                  "reused_completed_raw": report.get("reused_completed_raw", False),
                  "locks": inputs["locks"],
                  "source_environment": report["environment"],
                  "source_generation_report_sha256": digest(source / "generation.json"),
                  "source_swanlab_run_id": report["swanlab_run_id"],
                  "source_fasta_sha256": report["fasta_sha256"],
                  "score_cache": "All unique raw candidates; five scores and whole-reference audit frozen",
                  "note": "Model candidates; no activity or safety measurement"}
    save_json(output / "provenance.json", provenance)
    files = sorted(path for path in output.rglob("*") if path.is_file())
    manifest = {"format": "amp-frozen-pool-v1", "sampling_seed": inputs["seed"],
                "counts": complete["counts"], "raw_sha256": complete["sha256"],
                "sha256": {str(path.relative_to(output)): digest(path) for path in files}}
    save_json(output / "manifest.json", manifest)
    load_pool(output, verify_assets(output), load_reference(output / "reference"))
    archive.parent.mkdir(parents=True, exist_ok=True)
    with archive.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as handle:
            for path in sorted(path for path in output.rglob("*") if path.is_file()):
                info = handle.gettarinfo(str(path), arcname=str(path.relative_to(output)))
                info.uid = info.gid = info.mtime = 0
                info.uname = info.gname = ""
                info.mode = 0o644
                with path.open("rb") as stream:
                    handle.addfile(info, stream)
    descriptor = {"url": "https://github.com/nicetone9/amp-challenge-submission/releases/download/frozen-pool-seed42-20260930/" + archive.name,
                  "sha256": digest(archive), "manifest_sha256": digest(output / "manifest.json"),
                  "counts": complete["counts"], "raw_sha256": complete["sha256"]}
    save_json(archive.with_suffix(".json"), descriptor)
    print(json.dumps(descriptor), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--assets", type=Path, default=Path("checkpoint"))
    parser.add_argument("--output", type=Path, default=Path("frozen-pool"))
    parser.add_argument("--archive", type=Path, default=Path("work/release/frozen-pool-seed42.tar.gz"))
    args = parser.parse_args()
    freeze(args.source, args.assets, args.output, args.archive)
