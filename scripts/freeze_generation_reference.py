"""Freeze private reference-only inputs for fresh generation; never copy candidates."""
import argparse
import json
from pathlib import Path
import pandas as pd
from amp_submission.calibration import fingerprint
from amp_submission.hard_select import METRICS
from amp_submission.io import digest, read_fasta, verify_assets, write_fasta
from amp_submission.prepare_reference import save_json
from amp_submission.score_stage import load_scores


def freeze(work, assets, output):
    work, assets, output = map(Path, (work, assets, output))
    if output.exists():
        raise ValueError("Use a new immutable reference directory")
    verify_assets(assets)
    prepared = json.loads((work / "prepare-summary.json").read_text())
    for name in ("reference.csv", "training-sequences.csv"):
        if digest(work / name) != prepared["artifacts"][name]:
            raise ValueError("Reference/training input changed: " + name)
    reference, sources = load_scores(work, "reference")
    anchor = reference[reference.reference_role.eq("anchor")]
    calibration = reference[reference.reference_role.eq("calibration")]
    if set(anchor.cluster_id) & set(calibration.cluster_id) or set(anchor.sequence) & set(calibration.sequence):
        raise ValueError("Anchor/calibration leakage")
    motifs = json.loads((work / "motifs/frozen.json").read_text())
    if motifs["reference_sha256"] != digest(work / "reference.csv"):
        raise ValueError("Motifs use different reference")
    output.mkdir(parents=True)
    save_json(output / "reference-scores.json",
              {m: calibration[m].tolist() for m in METRICS})
    write_fasta(output / "anchors.fasta", anchor.sequence.tolist())
    known = sorted(set(reference.sequence) | set(read_fasta(assets / "antibacterial.fasta")))
    write_fasta(output / "known-amp.fasta", known)
    training = pd.read_csv(work / "training-sequences.csv").sequence
    (output / "training-hashes.txt").write_text("".join(h+"\n" for h in sorted({fingerprint(s) for s in training})))
    save_json(output / "motifs.json", motifs)
    save_json(output / "provenance.json",
              {"reference_sha256": digest(work / "reference.csv"),
               "training_sha256": digest(work / "training-sequences.csv"),
               "official_sha256": digest(assets / "antibacterial.fasta"),
               "motif_sha256": digest(work / "motifs/frozen.json"),
               "score_sources": sources, "anchors": len(anchor), "calibration": len(calibration),
               "known_amp_count": len(known), "candidate_data_included": False,
               "warning": "Private reference inputs, not authorized for public redistribution."})
    save_json(output / "manifest.json",
              {"format": 1, "sha256": {p.name: digest(p) for p in sorted(output.iterdir())}})
    verify_assets(output)
    print(json.dumps({"reference_manifest_sha256": digest(output / "manifest.json"),
                      "known_amp_count": len(known), "candidate_data_included": False}))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("work/six-metrics-600k"))
    p.add_argument("--assets", type=Path, default=Path("checkpoint"))
    p.add_argument("--output", type=Path, default=Path("checkpoint/generation-reference"))
    a = p.parse_args()
    freeze(a.work, a.assets, a.output)
