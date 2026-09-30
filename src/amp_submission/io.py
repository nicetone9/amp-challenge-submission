import hashlib
import json
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_fasta(path):
    sequences, parts = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if parts:
                sequences.append("".join(parts))
            parts = []
        elif line.strip():
            parts.append(line.strip())
    if parts:
        sequences.append("".join(parts))
    return sequences

def write_fasta(path, sequences):
    Path(path).write_text("".join(f">seq_{i:06d}\n{s}\n" for i, s in enumerate(sequences)))

def verify_assets(root):
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    for relative, expected in manifest["sha256"].items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file() or digest(path) != expected:
            raise ValueError(f"Missing or changed asset: {relative}")
    return manifest
