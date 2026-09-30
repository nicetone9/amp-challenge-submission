"""Official-reference hard novelty audit, separate from percentile qualification."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import fingerprint
from .generate import too_similar
from .io import digest, read_fasta
from .prepare_reference import save_json
from .score_stage import cache_scores


def audit_batch(sequences, references):
    return {"official_reference_pass": np.asarray(
        [not too_similar(sequence, references) for sequence in sequences], dtype=float)}


def run(work, assets):
    work, assets = Path(work), Path(assets)
    prepared = json.loads((work / "prepare-summary.json").read_text())
    candidate_path = work / "candidates.csv"
    if digest(candidate_path) != prepared["artifacts"]["candidates.csv"]:
        raise ValueError("Changed candidate table")
    reference_path = assets / "antibacterial.fasta"
    manifest = json.loads((assets / "manifest.json").read_text())
    expected = manifest["sha256"].get("antibacterial.fasta")
    if expected is None or digest(reference_path) != expected:
        raise ValueError("Missing or changed official reference")
    references = read_fasta(reference_path)
    if not references:
        raise ValueError("Empty official reference")
    identity = {"candidates_sha256": digest(candidate_path),
                "official_reference_sha256": expected,
                "code_sha256": digest(__file__),
                "comparison_code_sha256": digest(Path(__file__).with_name("generate.py")),
                "cache_code_sha256": digest(Path(__file__).with_name("score_stage.py")),
                "pixi_lock_sha256": digest("pixi.lock"),
                "max_ratio": .8, "batch_size": 512}
    target = work / "official-audit"
    state = target / "identity.json"
    if state.exists() and json.loads(state.read_text()) != identity:
        raise ValueError("Official audit fingerprint changed: use a new work directory")
    save_json(state, identity)
    candidates = pd.read_csv(candidate_path, keep_default_na=False)
    result = cache_scores(candidates.sequence.tolist(), target / "cache", identity,
                          lambda sequences: audit_batch(sequences, references))
    passed = result["official_reference_pass"].astype(bool)
    output = pd.DataFrame({"sequence": candidates.sequence,
                           "official_reference_pass": passed})
    output.to_csv(target / "candidates.csv", index=False)
    report = {"status": "HARD_NOVELTY_AUDITED_NOT_SIX_CATEGORY_QUALIFIED",
              "candidate_count": len(candidates), "official_reference_count": len(references),
              "official_reference_pass": int(passed.sum()),
              "official_reference_fail": int((~passed).sum()),
              "precheck_and_official_pass": int((passed & candidates.hard_precheck).sum()),
              "identity_sha256": fingerprint(identity),
              "result_sha256": digest(target / "candidates.csv"),
              "note": "Whole-library ratio <= 0.8 is the conservative project rule. "
                      "This does not replace AMP percentile or collection gates."}
    save_json(target / "complete.json", report)
    print(json.dumps(report), flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=Path("work/six-metrics"))
    parser.add_argument("--assets", type=Path, default=Path("checkpoint"))
    args = parser.parse_args()
    run(args.work, args.assets)


if __name__ == "__main__":
    main()
