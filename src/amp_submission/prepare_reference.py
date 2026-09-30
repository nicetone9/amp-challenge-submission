"""Audit the frozen milestone; never repurpose held-out sequences for selection."""
import argparse
import json
from collections import Counter
from pathlib import Path
import pandas as pd
from .calibration import fingerprint, registry_payload, split_clusters
from .io import AA, digest, read_fasta

def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)

def prepare(manifest_path, audit_path, candidate_root, output, seed=42, expanded_pool=None):
    output = Path(output)
    audit = json.loads(Path(audit_path).read_text())
    if digest(manifest_path) != audit["manifest_sha256"]:
        raise ValueError("Milestone manifest changed")
    for name, expected in audit["label_files"].items():
        if digest(name) != expected:
            raise ValueError("Step1 annotation input changed: " + name)
    paths = {arch: Path(candidate_root) / arch / "library.fasta" for arch in ("vq", "dima")}
    expanded = None
    if expanded_pool is not None:
        expanded_pool = Path(expanded_pool)
        expanded = json.loads((expanded_pool / "complete.json").read_text())
        pool_identity = json.loads((expanded_pool / "identity.json").read_text())
        if expanded["identity_sha256"] != fingerprint(pool_identity):
            raise ValueError("Expansion identity mismatch")
        if pool_identity["seed"] != seed:
            raise ValueError("Expansion/reference seed mismatch")
        if any(pool_identity["initial"][arch]["sha256"] != digest(paths[arch]) for arch in paths):
            raise ValueError("Expansion used different initial libraries")
    identity = {"manifest_sha256": audit["manifest_sha256"], "label_files": audit["label_files"],
                "seed": seed, "input_fasta_sha256": {k: digest(p) for k, p in paths.items()},
                "registry_sha256": fingerprint(registry_payload()), "code_sha256": digest(__file__),
                "expansion_sha256": digest(expanded_pool / "complete.json") if expanded else None}
    if (output / "identity.json").exists():
        if json.loads((output / "identity.json").read_text()) != identity:
            raise ValueError("Refusing changed inputs in frozen reference directory")
    frame = pd.read_csv(manifest_path, keep_default_na=False)
    if not frame.sequence.is_unique or frame.groupby("cluster_id").split.nunique().max() != 1:
        raise ValueError("Duplicate sequences or cross-split cluster leakage")
    if frame.cluster_id.eq("").any():
        raise ValueError("Missing cluster identifiers")
    train = frame[frame.split.eq("train")]
    core = train[train.stratum.eq("core")].copy()
    eligible = core.length.between(8, 40) & core.sequence.str.fullmatch("[" + AA + "]+")
    excluded = core.loc[~eligible].copy()
    core = core.loc[eligible].sort_values("sequence_key", kind="stable").copy()
    core["reference_role"] = split_clusters(core.cluster_id.tolist(), seed)
    core["chemical_state"] = "unknown_in_sequence_only_milestone"
    core["evidence"] = "frozen ZYL AMP train-positive or OmegAMP curated_amp"
    known = set(train.sequence)
    merged, seen = [], {}
    counts = Counter()
    for arch, path in paths.items():
        sequences = read_fasta(path)
        if len(sequences) != 50000 or len(set(sequences)) != 50000:
            raise ValueError("Initial branch must contain exactly 50000 unique sequences: " + arch)
        if expanded:
            needed = expanded["additional_counts"][arch]
            extra = []
            index = 0
            while len(extra) < needed:
                batch = json.loads((expanded_pool / arch / f"batch-{index:07d}.json").read_text())
                expected = fingerprint({"identity": pool_identity, "arch": arch, "index": index,
                                        "batch_size": pool_identity["batch_size"]})
                if (batch["identity"] != expected or
                        batch["sequences_sha256"] != fingerprint(batch["sequences"]) or
                        len(batch["sequences"]) != pool_identity["batch_size"]):
                    raise ValueError("Corrupt expansion batch")
                extra.extend(batch["sequences"])
                index += 1
            sequences.extend(extra[:needed])
        counts[arch + "_raw"] = len(sequences)
        for sequence in sequences:
            reasons = []
            if not 8 <= len(sequence) <= 40 or not set(sequence) <= set(AA):
                reasons.append("unsupported_sequence")
            if sequence and max(Counter(sequence).values()) / len(sequence) > .6:
                reasons.append("single_residue_fraction_gt_0.6")
            if sequence in known:
                reasons.append("exact_training_hit")
            if sequence in seen:
                sources = set(seen[sequence]["source"].split("|"))
                if arch not in sources:
                    counts["cross_model_duplicates"] += 1
                else:
                    counts["within_model_duplicates"] += 1
                seen[sequence]["source"] = "|".join(sorted(sources | {arch}))
                seen[sequence]["raw_occurrences"] += 1
                continue
            row = {"sequence": sequence, "source": arch, "sequence_sha256": fingerprint(sequence),
                   "hard_precheck": not reasons, "failure_reasons": ";".join(reasons), "raw_occurrences": 1}
            seen[sequence] = row
            merged.append(row)
    output.mkdir(parents=True, exist_ok=True)
    core.to_csv(output / "reference.csv", index=False)
    excluded.to_csv(output / "excluded-reference.csv", index=False)
    pd.DataFrame(merged).to_csv(output / "candidates.csv", index=False)
    # Private asset used only for exact novelty, never uploaded with a report.
    pd.DataFrame({"sequence": sorted(known)}).to_csv(output / "training-sequences.csv", index=False)
    save_json(output / "identity.json", identity)
    summary = {"status": "PREPARED_NOT_CALIBRATED", **identity,
               "reference_train_core_all_lengths": len(core) + len(excluded),
               "reference_eligible": len(core), "excluded_length_or_alphabet": len(excluded),
               "reference_roles": core.reference_role.value_counts().to_dict(),
               "reference_clusters": core.groupby("reference_role").cluster_id.nunique().to_dict(),
               "pool_raw": dict(counts), "pool_unique": len(merged),
               "pool_precheck_pass": sum(row["hard_precheck"] for row in merged),
               "hard_precheck_failures": dict(Counter(reason for row in merged
                   for reason in row["failure_reasons"].split(";") if reason)),
               "data_boundary": "Only frozen milestone train-core. Not all newer Step1 databases.",
               "limitations": ["Chemical states absent from milestone; no modification inference from letters",
                              "Additional Step1 positives outside frozen milestone require an audited manifest",
                              "Official-reference similarity audit runs at scoring/selection stage"],
               "artifacts": {name: digest(output / name) for name in
                             ("reference.csv", "excluded-reference.csv", "candidates.csv", "training-sequences.csv")}}
    save_json(output / "prepare-summary.json", summary)
    save_json(output / "metric-registry.json", registry_payload())
    print(json.dumps({k: v for k, v in summary.items() if k not in ("label_files", "limitations")}))
    return summary

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--audit", type=Path, required=True)
    p.add_argument("--candidate-root", type=Path, default=Path("candidates"))
    p.add_argument("--output", type=Path, default=Path("work/six-metrics"))
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--expanded-pool", type=Path)
    a = p.parse_args()
    prepare(a.manifest, a.audit, a.candidate_root, a.output, a.seed, a.expanded_pool)

if __name__ == "__main__":
    main()
