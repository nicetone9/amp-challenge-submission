"""Competition entry point. All arguments have defaults; missing assets fail closed."""
import argparse
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Indel
from .io import AA, digest, read_fasta, write_fasta, verify_assets

def too_similar(sequence, references):
    return process.extractOne(sequence, references, scorer=Indel.normalized_similarity,
                              score_cutoff=np.nextafter(.8, 1.0)) is not None

def rank_top(pool, joint, count=100):
    selected, available, near = [], np.ones(len(pool), dtype=bool), np.zeros(len(pool))
    for _ in range(count):
        utility = np.asarray(joint) - .15 * near
        utility[~available] = -np.inf
        index = int(np.argmax(utility))
        selected.append(pool[index])
        available[index] = False
        near = np.maximum(near, [Indel.normalized_similarity(s, pool[index]) for s in pool])
    return selected

def validate(library, top, references, library_size=50000, top_size=100):
    if len(library) != library_size or len(set(library)) != library_size:
        raise ValueError("Library count or uniqueness check failed")
    if len(top) != top_size or len(set(top)) != top_size or not set(top) <= set(library):
        raise ValueError("Top count, uniqueness or membership check failed")
    if any(not 8 <= len(s) <= 50 or not set(s) <= set(AA) for s in library):
        raise ValueError("Noncanonical sequence or invalid length")
    if any(too_similar(s, references) for s in library):
        raise ValueError("Official-reference similarity exceeds 0.8")

def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--arch", choices=("vq", "dima"), default="vq")
    p.add_argument("--assets", type=Path, default=Path("checkpoint"))
    p.add_argument("--oracles", type=Path, default=Path("external"))
    p.add_argument("--output", type=Path, default=Path("generate"))
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    p.add_argument("--threads", type=int, default=2)
    return p

def main():
    args = parser().parse_args()
    # Verify every frozen asset before any generation or output mutation.
    manifest = verify_assets(args.assets)
    oracle_manifest = verify_assets(args.oracles)
    import torch
    from .sampling import OBJECTIVES, seed_all, load_models, sample
    from .scoring import predict_external, rewards
    torch.set_num_threads(args.threads)
    seed_all(args.seed)
    device = ("cuda" if torch.cuda.is_available() else "cpu") if args.device == "auto" else args.device
    policies, reference, bundle = load_models(args.assets, args.arch, device)
    pdf = np.asarray(json.loads((args.assets / "length-pdf.json").read_text()), dtype=float)
    pdf[41:] = 0
    pdf /= pdf.sum()
    references = read_fasta(args.assets / "antibacterial.fasta")
    calibration = json.loads((args.assets / "calibration.json").read_text())["scores"]
    # Model construction consumes RNG; match the evaluated reseed before sampling.
    seed_all(args.seed)
    library, seen, counts = [], set(), Counter()
    for batch_id in range((500000 + 127) // 128):
        objective = OBJECTIVES[batch_id % len(OBJECTIVES)]
        n = min(128, 500000 - counts["raw"])
        sequences = sample(args.arch, policies[objective], reference, bundle,
                           np.random.choice(len(pdf), n, p=pdf).tolist(), device)
        for sequence in sequences:
            counts["raw"] += 1
            if len(library) == 50000:
                counts["unused_after_target"] += 1
                continue
            if sequence in seen:
                counts["duplicate"] += 1
                continue
            seen.add(sequence)
            if not 8 <= len(sequence) <= 40 or not set(sequence) <= set(AA):
                counts["invalid"] += 1
                continue
            if max(Counter(sequence).values()) / len(sequence) > .6:
                counts["low_complexity"] += 1
                continue
            if too_similar(sequence, references):
                counts["reference_similar"] += 1
                continue
            library.append(sequence)
        if len(library) == 50000:
            break
    if len(library) != 50000:
        raise RuntimeError(f"Candidate budget exhausted: {len(library)}/50000; {dict(counts)}")
    indices = np.random.default_rng(4202048).choice(50000, 2048, replace=False)
    pool = [library[int(i)] for i in indices]
    scores = predict_external(pool, args.oracles, device)
    top = rank_top(pool, rewards(scores, calibration)["joint"])
    validate(library, top, references)
    args.output.mkdir(parents=True, exist_ok=True)
    # Repeated invocation in the same directory is supported by the official runner.
    # Failed generation leaves previous FASTAs intact; each final file is atomically replaced.
    with tempfile.TemporaryDirectory(prefix=".amp-staging-", dir=args.output.parent) as temp:
        staging = Path(temp)
        write_fasta(staging / "library.fasta", library)
        write_fasta(staging / "top.fasta", top)
        report = {"arch": args.arch, "seed": args.seed, "device": device,
                  "counts": dict(counts), "assets": digest(args.assets / "manifest.json"),
                  "oracles": digest(args.oracles / "manifest.json"),
                  "library_sha256": digest(staging / "library.fasta"),
                  "top_sha256": digest(staging / "top.fasta"),
                  "status": "component_checks_passed_not_organizer_acceptance"}
        (staging / "generation.json").write_text(json.dumps(report, indent=2) + "\n")
        for name in ("library.fasta", "top.fasta", "generation.json"):
            os.replace(staging / name, args.output / name)
    print(json.dumps(report))

if __name__ == "__main__":
    main()
