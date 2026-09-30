"""Audit the length-hard-only fixed-pool delivery and its independent repeat."""
import json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from amp_submission.calibration import fingerprint
from amp_submission.generate import validate
from amp_submission.hard_select import METRICS
from amp_submission.io import digest, read_fasta
from amp_submission.motif_select import motif_top
from amp_submission.prepare_reference import save_json
from amp_submission.score_stage import load_scores

ROOT = Path("work/length-hard-only-seed42")
REPEAT = Path("work/length-hard-only-seed42-repeat")
BASE = Path("work/length-hard-only-seed42-base")
WORK = Path("work/six-metrics-600k")


def cdf(reference, values):
    reference = np.sort(np.asarray(reference, float))
    values = np.asarray(values, float)
    return (np.searchsorted(reference, values, side="left") +
            np.searchsorted(reference, values, side="right")) / (2 * len(reference))


def calibrated(reference, values, typical):
    if not typical:
        return cdf(reference, values)
    f_ref, f_values = cdf(reference, reference), cdf(reference, values)
    return cdf(2 * np.minimum(f_ref, 1 - f_ref), 2 * np.minimum(f_values, 1 - f_values))


def length_summary(sequences):
    lengths = np.asarray([len(s) for s in sequences])
    return {"count": len(sequences), "min": int(lengths.min()), "max": int(lengths.max()),
            "mean": float(lengths.mean()), "median": float(np.median(lengths)),
            "counts": dict(sorted(Counter(map(int, lengths)).items())),
            "bins": {f"{lo}-{hi}": int(((lengths >= lo) & (lengths <= hi)).sum())
                     for lo, hi in ((8, 15), (16, 20), (21, 25), (26, 30), (31, 40), (41, 50))}}


def main():
    base = json.loads((BASE / "complete.json").read_text())
    complete = json.loads((ROOT / "complete.json").read_text())
    assert base["quality_metrics"] == list(METRICS)
    assert base["descriptive_only"] == ["length", "molecular_weight"]
    assert not {"length", "molecular_weight"} & set(METRICS)
    assert digest(WORK / "candidates.csv") == base["identity"]["candidate_sha256"]
    assert digest(WORK / "reference.csv") == base["identity"]["reference_sha256"]
    sequences, tables, hashes = {}, {}, {}
    for name, count in (("library", 50000), ("top", 100)):
        path = ROOT / (name + ".fasta")
        assert path.read_bytes() == (REPEAT / path.name).read_bytes()
        assert digest(path) == complete["fasta_sha256"][name]
        hashes[name] = digest(path)
        sequences[name] = read_fasta(path)
        assert len(sequences[name]) == len(set(sequences[name])) == count
        tables[name] = pd.read_csv(ROOT / (name + "-scores.csv"), keep_default_na=False)
        assert tables[name].sequence.tolist() == sequences[name]
        assert Counter(tables[name].tier) == {"high": count // 2, "low": count // 2}
        assert not {"length_percentile", "molecular_weight_percentile"} & set(tables[name])
    assert set(sequences["top"]) <= set(sequences["library"])
    reference, _ = load_scores(WORK, "reference")
    candidates, _ = load_scores(WORK, "candidates")
    reference = reference[reference.reference_role.eq("calibration")]
    selected = candidates.set_index("sequence").loc[sequences["library"]]
    assert selected.hard_precheck.all()
    values = np.column_stack([calibrated(reference[m], selected[m], m != "anchor_novelty")
                              for m in METRICS])
    assert np.isfinite(values).all() and (values >= base["selected_threshold"]).all()
    np.testing.assert_allclose(tables["library"][[m + "_percentile" for m in METRICS]],
                               values, rtol=0, atol=1e-12)
    ranking = (values[:, :-1].mean(axis=1) + values[:, -1]) / 2
    np.testing.assert_allclose(tables["library"].ranking, ranking, rtol=0, atol=1e-12)
    np.testing.assert_allclose(tables["library"].weakest, values.min(axis=1), rtol=0, atol=1e-12)
    motifs = json.loads((WORK / "motifs/frozen.json").read_text())["motifs"]
    short = sorted(m for m in motifs if 3 <= len(m) <= 4)
    for table in tables.values():
        for row in table.to_dict("records"):
            assert row["short_motif_support"] == sum(m in row["sequence"] for m in short)
            assert row["motif_support"] == sum(m in row["sequence"] for m in motifs)
    top_rows = tables["top"].to_dict("records")
    rank_key = lambda r: (-r["ranking"], -r["weakest"], -r["motif_support"], fingerprint(r["sequence"]))
    assert [r["sequence"] for r in sorted(top_rows, key=rank_key)] == sequences["top"]
    assert [r["sequence"] for r in motif_top(tables["library"].to_dict("records"), 100, 10)] == sequences["top"]
    quota = tables["top"][tables["top"].selection_reason.eq("short_motif_quota")]
    assert len(quota) == 10 and (quota.short_motif_support > 0).all()
    assert Counter(quota.tier) == {"high": 5, "low": 5}
    known = sorted(set(read_fasta("checkpoint/antibacterial.fasta")) |
                   set(pd.read_csv(WORK / "reference.csv").sequence))
    print("Independent full-library known-reference hard audit starting", flush=True)
    validate(sequences["library"], sequences["top"], known, 50000, 100)
    report = {"status": "LENGTH_HARD_ONLY_FIXED_POOL_DOUBLE_RUN_VERIFIED",
              "seed": 42, "scope": "Frozen-pool selection only; not two fresh model generations",
              "quality_metrics": list(METRICS), "descriptive_only": ["length", "molecular_weight"],
              "hard_length_bounds": [8, 50], "selected_threshold": base["selected_threshold"],
              "threshold_curve": base["threshold_curve"], "fasta_sha256": hashes,
              "byte_identical": True, "independent_percentile_and_rank_recalculation": True,
              "whole_library_known_reference_ratio_le_0_8": True, "known_reference_count": len(known),
              "all_six_categories_passed": False,
              "motif_length_range": [3, 4], "motif_quota": 10,
              "statistics": {name: {**length_summary(sequences[name]),
                  "source_counts": dict(Counter(tables[name].source)),
                  "tier_counts": dict(Counter(tables[name].tier)),
                  "clusters": int(tables[name].cluster.nunique()),
                  "motif_supported": int((tables[name].short_motif_support > 0).sum())}
                  for name in sequences},
              "inputs": base["identity"],
              "audit_code_sha256": digest(__file__),
              "warning": "Internal five-metric ranking; no official aggregation or experimental efficacy claim."}
    save_json("reports/length-hard-only-seed42-verification.json", report)
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
