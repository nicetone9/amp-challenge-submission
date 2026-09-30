"""Select a hard-compliant library using an explicitly partial quality panel."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import REGISTRY, fingerprint, percentile
from .io import AA, digest, read_fasta, write_fasta
from .prepare_reference import save_json
from .score_stage import load_scores
from .selection import assign_tiers, cluster_sequences, select

METRICS = ("length", "charge_ph7", "gravy", "aromaticity",
           "isoelectric_point", "molecular_weight", "anchor_novelty")
THRESHOLDS = (.50, .45, .40, .35, .30, .25)


def sequence_hard_mask(sequences, official_reference):
    return np.asarray([8 <= len(s) <= 50 and set(s) <= set(AA)
                       and s not in official_reference for s in sequences])


def quality_arrays(reference, candidates):
    """Equal weight within categories and between the two available categories."""
    registered = {m.name: m for m in REGISTRY}
    values = np.column_stack([percentile(reference[n], candidates[n],
                                        registered[n].direction) for n in METRICS])
    groups = {}
    for i, name in enumerate(METRICS):
        groups.setdefault(registered[name].category, []).append(i)
    ranking = np.mean([values[:, indexes].mean(axis=1) for indexes in groups.values()], axis=0)
    return values, ranking, values.min(axis=1)


def threshold_curve(values, hard_mask):
    return [{"threshold": t, "eligible": int((hard_mask & np.isfinite(values).all(axis=1)
                                             & (values >= t).all(axis=1)).sum())}
            for t in THRESHOLDS]


def run(work, output, assets, count=50000, top_k=100, seed=42, export=False):
    work, output, assets = Path(work), Path(output), Path(assets)
    if count < 1 or not 1 <= top_k <= count:
        raise ValueError("Require 1 <= top-k <= n-sequences")
    if (output / "complete.json").exists():
        raise ValueError("Do not overwrite completed output; choose a new output directory")
    prepared = json.loads((work / "prepare-summary.json").read_text())
    for name, expected in prepared["artifacts"].items():
        if digest(work / name) != expected:
            raise ValueError("Prepared artifact changed: " + name)
    reference, ref_sources = load_scores(work, "reference")
    candidates, sources = load_scores(work, "candidates")
    reference = reference[reference.reference_role.eq("calibration")]
    audit = json.loads((work / "official-audit/complete.json").read_text())
    if digest(work / "official-audit/candidates.csv") != audit["result_sha256"]:
        raise ValueError("Changed official audit results")
    audit_identity = json.loads((work / "official-audit/identity.json").read_text())
    if audit_identity["candidates_sha256"] != digest(work / "candidates.csv"):
        raise ValueError("Official audit belongs to a different candidate pool")
    if digest(assets / "antibacterial.fasta") != audit_identity["official_reference_sha256"]:
        raise ValueError("Official reference changed")
    flags = pd.read_csv(work / "official-audit/candidates.csv")
    if len(flags) != len(candidates) or set(flags.sequence) != set(candidates.sequence):
        raise ValueError("Incomplete official audit coverage")
    candidates = candidates.merge(flags, on="sequence", validate="one_to_one")
    if not candidates.sequence.is_unique:
        raise ValueError("Duplicate candidates")
    sequences = candidates.sequence.tolist()
    hard = sequence_hard_mask(sequences, set(read_fasta(assets / "antibacterial.fasta")))
    base = hard & candidates.hard_precheck.to_numpy(bool) & candidates.official_reference_pass.to_numpy(bool)
    values, ranking, weakest = quality_arrays(reference, candidates)
    lengths = candidates.sequence.str.len()
    report = {
        "status": "AUDIT_ONLY", "seed": seed, "candidate_count": len(candidates),
        "sequence_hard_pass": int(hard.sum()),
        "standard_amino_acids_pass": int(candidates.sequence.map(lambda s: set(s) <= set(AA)).sum()),
        "length_8_50_pass": int(lengths.between(8, 50).sum()),
        "length_min": int(lengths.min()), "length_max": int(lengths.max()),
        "length_bins": {f"{lo}-{hi}": int(lengths.between(lo, hi).sum())
                        for lo, hi in ((8, 15), (16, 25), (26, 40), (41, 50))},
        "conservative_hard_and_precheck_pass": int(base.sum()),
        "quality_metrics": list(METRICS), "threshold_curve": threshold_curve(values, base),
        "quality_complete": False, "all_six_categories_passed": False,
        "chemical_design": "Linear, unmodified, free N and C termini; not an experimental measurement.",
        "pending": ["Full potency and safety prediction", "ProtT5 and collection metrics",
                    "Validated peptide synthesizability", "Full generation reproducibility and package validation"],
        "identity": {"candidate_sha256": digest(work / "candidates.csv"),
                     "reference_sha256": digest(work / "reference.csv"),
                     "reference_score_sources": ref_sources, "candidate_score_sources": sources,
                     "official_result_sha256": audit["result_sha256"],
                     "code_sha256": digest(__file__),
                     "selection_code_sha256": digest(Path(__file__).with_name("selection.py")),
                     "pixi_lock_sha256": digest("pixi.lock")},
    }
    save_json(output / "audit.json", report)
    print(json.dumps(report), flush=True)
    if not export:
        return report
    motif = json.loads((work / "motifs/frozen.json").read_text())
    if motif["reference_sha256"] != digest(work / "reference.csv"):
        raise ValueError("Motif reference changed")
    report["identity"]["motif_sha256"] = digest(work / "motifs/frozen.json")
    attempts = []
    chosen = None
    for item in report["threshold_curve"]:
        t = item["threshold"]
        if item["eligible"] < count:
            continue
        mask = base & np.isfinite(values).all(axis=1) & (values >= t).all(axis=1)
        indices = np.flatnonzero(mask)
        seqs = candidates.iloc[indices].sequence.tolist()
        clusters, command = cluster_sequences(seqs, output / f"cluster-{t:.2f}")
        rows = [{"sequence": sequences[i], "source": candidates.iloc[i].source,
                 "cluster": clusters[sequences[i]], "passed": True,
                 "ranking": float(ranking[i]), "weakest": float(weakest[i]),
                 "motif_support": sum(m in sequences[i] for m in motif["motifs"]),
                 **{name + "_percentile": float(values[i, j]) for j, name in enumerate(METRICS)}}
                for i in indices]
        tiers = assign_tiers(rows, seed)
        capacities = {tier: sum(r["tier"] == tier for r in tiers) for tier in ("high", "low")}
        attempts.append({"threshold": t, "tier_capacity": capacities, "command": command})
        if capacities["high"] < count - count // 2 or capacities["low"] < count // 2:
            continue
        chosen = select(tiers, count)
        top = select(chosen, top_k)
        report["selected_threshold"] = t
        report["cluster_count"] = len(set(clusters.values()))
        break
    if chosen is None:
        report.update(status="INSUFFICIENT_WITHIN_APPROVED_THRESHOLD_GRID", selection_attempts=attempts)
        save_json(output / "audit.json", report)
        raise ValueError("Cannot fill quotas at thresholds >= 0.25; no automatic relaxation")
    library_sequences = [r["sequence"] for r in chosen]
    top_sequences = [r["sequence"] for r in top]
    assert len(set(library_sequences)) == count
    assert len(set(top_sequences)) == top_k and set(top_sequences) <= set(library_sequences)
    output.mkdir(parents=True, exist_ok=True)
    # Stage outputs and promote only after all in-memory gates have passed.
    for name, seqs in (("library", library_sequences), ("top", top_sequences)):
        staged = output / (name + ".fasta.partial")
        write_fasta(staged, seqs)
        staged.replace(output / (name + ".fasta"))
    pd.DataFrame(chosen).to_csv(output / "library-scores.csv", index=False)
    pd.DataFrame(top).to_csv(output / "top-scores.csv", index=False)
    report.update(status="HARD_SEQUENCE_COMPLIANT_PARTIAL_QUALITY_LIBRARY",
                  n_sequences=count, top_k=top_k, selection_attempts=attempts,
                  source_counts=pd.Series([r["source"] for r in chosen]).value_counts().to_dict(),
                  tier_counts=pd.Series([r["tier"] for r in chosen]).value_counts().to_dict(),
                  fasta_sha256={name: digest(output / (name + ".fasta")) for name in ("library", "top")},
                  note="Internal ranking uses only the seven named metrics, not an official score. "
                       "The passed flag means this partial panel and hard gates only, not all six categories.")
    save_json(output / "complete.json", report)
    print(json.dumps(report), flush=True)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("work/six-metrics-600k"))
    p.add_argument("--output", type=Path, default=Path("work/hard-library-seed42"))
    p.add_argument("--assets", type=Path, default=Path("checkpoint"))
    p.add_argument("--n-sequences", type=int, default=50000)
    p.add_argument("--top-k", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--export", action="store_true")
    a = p.parse_args()
    run(a.work, a.output, a.assets, a.n_sequences, a.top_k, a.seed, a.export)


if __name__ == "__main__":
    main()
