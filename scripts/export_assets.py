"""Export trusted local experiment checkpoints. Never run this on untrusted pickle files."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import torch

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("experiment", type=Path)
    p.add_argument("--output", type=Path, default=Path("checkpoint"))
    a = p.parse_args()
    source, out = a.experiment.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    provenance = {}
    def export(relative, checkpoint, expected=None):
        if expected and digest(checkpoint) != expected:
            raise ValueError(f"Source checkpoint changed: {checkpoint}")
        target = out / relative
        if target.exists():
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        state = torch.load(checkpoint, map_location="cpu", weights_only=False)
        torch.save({k: v.detach().cpu() for k, v in state["model"].items()}, target)
        provenance[relative] = {"source_sha256": digest(checkpoint),
                               "source_relative": str(checkpoint.relative_to(source))}
    for arch in ("vq", "dima"):
        config = json.loads((source / f"challenge_20260929/deliverables/{arch}/config.json").read_text())
        for objective, row in config["selection"].items():
            original = Path(row["path"])
            # Paths are confined to the explicitly provided trusted experiment root.
            if not original.resolve().is_relative_to(source):
                raise ValueError("Checkpoint is outside trusted experiment root")
            export(f"{arch}/{objective}.pt", original, row["sha256"])
            provenance[f"{arch}/{objective}.pt"].update(
                {k: row[k] for k in ("run_id", "seed", "selected_step", "validation_score")})
    for relative, origin in [
        ("vq/codec.pt", "checkpoints/vq-codec/best.pt"),
        ("dima/decoder.pt", "checkpoints/dima-decoder/best.pt"),
        ("dima/reference.pt", "extension_20260928/dima/best.pt")]:
        cp = source / origin
        completion = json.loads((cp.parent / "complete.json").read_text())
        export(relative, cp, completion["best_sha256"])
    for name, origin in [
        ("esm-config.json", "models/esm2_650m/config.json"),
        ("normalization.json", "cache_esm/normalization.json"),
        ("calibration.json", "challenge_20260929/calibration.json"),
        ("antibacterial.fasta", "challenge_20260929/official/data/antibacterial.fasta")]:
        shutil.copyfile(source / origin, out / name)
    audit = json.loads((source / "prepared/audit.json").read_text())
    (out / "length-pdf.json").write_text(json.dumps(audit["weighted_train_length_pdf"]) + "\n")
    hashes = {str(f.relative_to(out)): digest(f) for f in sorted(out.rglob("*"))
              if f.is_file() and f.name not in ("manifest.json", "README.md")}
    (out / "manifest.json").write_text(json.dumps(
        {"format": 1, "provenance": provenance, "sha256": hashes}, indent=2) + "\n")
    print(json.dumps({"exported_files": len(hashes), "manifest_sha256": digest(out / "manifest.json")}))

if __name__ == "__main__":
    main()
