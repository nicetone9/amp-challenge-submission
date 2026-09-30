"""Rerank an unchanged hard-compliant library with a short-motif quota."""
import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from .generate import too_similar
from .hard_select import METRICS
from .io import digest, read_fasta, write_fasta
from .motifs import scan
from .prepare_reference import save_json
from .selection import rank_key, select


def motif_top(rows, count=100, quota=10):
    if not 0 <= quota <= count or quota % 2 or count % 2:
        raise ValueError("Require even top count and even motif quota in range")
    if len({r["sequence"] for r in rows}) != len(rows):
        raise ValueError("Duplicate library rows")
    anchors = select([r for r in rows if r["short_motif_support"]], quota) if quota else []
    selected = {r["sequence"] for r in anchors}
    remainder = select([r for r in rows if r["sequence"] not in selected], count-quota) if count>quota else []
    result = [{**r, "selection_reason": "short_motif_quota"} for r in anchors]
    result += [{**r, "selection_reason": "cluster_quality_quota"} for r in remainder]
    return sorted(result, key=rank_key)


def step1_support(root, allowed, motifs):
    """Use only frozen-reference members also safe in current Step 1 train."""
    root = Path(root)
    hard_path = root / "audit/training/global_hard_isolated_train.csv"
    hard = pd.read_csv(hard_path, low_memory=False)
    positive = set(hard.loc[hard.task.isin(["antimicrobial", "antibacterial"])
                           & hard.target.eq(1), "sequence"]) & allowed
    safe_sequences, hashes = set(), {str(hard_path): digest(hard_path)}
    for task in ("antimicrobial", "antibacterial"):
        path = root / "amp" / task / "train.csv"
        frame = pd.read_csv(path, keep_default_na=False, low_memory=False)
        hashes[str(path)] = digest(path)
        safe = frame.target.eq(1) & frame.sequence.isin(positive)
        for column in ("qmap_holdout", "fixed_test", "escape_test", "multipep_test",
                       "external_benchmark_test", "q_due_to_motif", "q_due_to_h30",
                       "q_due_to_cross_task", "paper_label_conflict"):
            if column in frame:
                safe &= ~frame[column].astype(str).str.lower().isin(["true", "1"])
        for column in ("original_assigned_split", "original_split"):
            if column in frame:
                safe &= ~frame[column].astype(str).str.lower().str.contains("test|quarantine|valid")
        safe_sequences.update(frame.loc[safe, "sequence"])
    return {"eligible_unique": len(safe_sequences), "input_sha256": hashes,
            "motif_support": {m: sum(m in s for s in safe_sequences) for m in motifs},
            "scope": "Intersection with frozen AMP training reference; no new reference expansion."}


def run(library, work, assets, step1, output, count=100, quota=10):
    library, work, assets, output = map(Path, (library, work, assets, output))
    if output.exists():
        raise ValueError("Use a new output directory; never overwrite a frozen delivery")
    original = json.loads((library / "complete.json").read_text())
    if digest(library / "library.fasta") != original["fasta_sha256"]["library"]:
        raise ValueError("Changed input library")
    frozen_path = work / "motifs/frozen.json"
    if digest(frozen_path) != original["identity"]["motif_sha256"]:
        raise ValueError("Changed motif reference")
    if digest(work / "reference.csv") != original["identity"]["reference_sha256"]:
        raise ValueError("Changed AMP reference")
    frozen = json.loads(frozen_path.read_text())
    motifs = sorted(m for m in frozen["motifs"] if 3 <= len(m) <= 4)
    ref = pd.read_csv(work / "reference.csv", keep_default_na=False)
    step1_evidence = step1_support(step1, set(ref.sequence), motifs)
    official_path = assets / "antibacterial.fasta"
    audit_identity = json.loads((work / "official-audit/identity.json").read_text())
    if digest(official_path) != audit_identity["official_reference_sha256"]:
        raise ValueError("Changed official reference")
    audit = json.loads((work / "official-audit/complete.json").read_text())
    if digest(work / "official-audit/candidates.csv") != audit["result_sha256"]:
        raise ValueError("Changed official audit")
    flags = pd.read_csv(work / "official-audit/candidates.csv").set_index("sequence")
    sequences = read_fasta(library / "library.fasta")
    table = pd.read_csv(library / "library-scores.csv")
    if table.sequence.tolist() != sequences or len(set(sequences)) != 50000:
        raise ValueError("Library and scores mismatch")
    values = table[[m+"_percentile" for m in METRICS]].to_numpy()
    if not np.isfinite(values).all() or not (values >= original["selected_threshold"]).all():
        raise ValueError("Input quality gates failed")
    if not flags.loc[sequences, "official_reference_pass"].all():
        raise ValueError("Official reference gate failed")
    extra = sorted(set(ref.sequence)-set(read_fasta(official_path)))
    failed = [s for s in sequences if too_similar(s, extra)]
    if failed:
        raise ValueError(f"{len(failed)} library sequences fail additional AMP reference gate; refill required")
    rows = table.to_dict("records")
    for row in rows:
        hits = scan(row["sequence"], motifs)
        row["short_motif_support"] = len({h["motif"] for h in hits})
        row["short_motifs"] = ";".join(sorted({h["motif"] for h in hits}))
        row["motif_positions_json"] = json.dumps(hits, sort_keys=True, separators=(",", ":"))
    top = motif_top(rows, count, quota)
    assert len(top)==count and len({r["sequence"] for r in top})==count
    assert sum(bool(r["short_motif_support"]) for r in top)>=quota
    assert sum(r["tier"]=="high" for r in top)==count//2
    output.mkdir(parents=True)
    shutil.copyfile(library / "library.fasta", output / "library.fasta")
    write_fasta(output / "top.fasta", [r["sequence"] for r in top])
    pd.DataFrame(rows).to_csv(output / "library-scores.csv", index=False)
    pd.DataFrame(top).to_csv(output / "top-scores.csv", index=False)
    report = {"status": "HARD_COMPLIANT_SHORT_MOTIF_QUOTA_PARTIAL_QUALITY",
              "all_six_categories_passed": False, "seed": original["seed"],
              "n_sequences": len(sequences), "top_k": count, "motif_quota": quota,
              "motif_length_range": [3,4], "reference_percentile": original["selected_threshold"],
              "library_unchanged": digest(output/"library.fasta")==original["fasta_sha256"]["library"],
              "motif_supported_library": sum(bool(r["short_motif_support"]) for r in rows),
              "motif_supported_top": sum(bool(r["short_motif_support"]) for r in top),
              "top_source_counts": pd.Series([r["source"] for r in top]).value_counts().to_dict(),
              "top_tier_counts": pd.Series([r["tier"] for r in top]).value_counts().to_dict(),
              "top_clusters": len({r["cluster"] for r in top}),
              "additional_amp_reference_count": len(extra), "additional_reference_fail": len(failed),
              "step1_evidence": step1_evidence,
              "database_evidence": {m: {source: int((ref.sequence.str.contains(m,regex=False) &
                  ref.source_membership.str.contains(source,regex=False)).sum())
                  for source in ("ZYL_dataset","OmegAMP")} for m in motifs},
              "input_sha256": {"library": digest(library/"library.fasta"),
                  "scores": digest(library/"library-scores.csv"), "motifs": digest(frozen_path),
                  "reference": digest(work/"reference.csv"), "official_reference": digest(official_path),
                  "code": digest(__file__), "selection_code": digest(Path(__file__).with_name("selection.py")),
                  "lock": digest("pixi.lock")},
              "fasta_sha256": {name: digest(output/(name+".fasta")) for name in ("library","top")},
              "warning": "Short motifs are overlapping support, not proof of activity. Full predictor, embedding, synthesis and fresh-generation checks remain incomplete."}
    save_json(output/"complete.json", report)
    print(json.dumps(report), flush=True)
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library",type=Path,default=Path("work/hard-library-seed42"))
    parser.add_argument("--work",type=Path,default=Path("work/six-metrics-600k"))
    parser.add_argument("--assets",type=Path,default=Path("checkpoint"))
    parser.add_argument("--step1",type=Path,default=Path("/home/peiranj/AMP/step1_retrain_external_20260917/data/ZYL_dataset"))
    parser.add_argument("--output",type=Path,default=Path("work/motif-library-seed42"))
    parser.add_argument("--top-k",type=int,default=100)
    parser.add_argument("--motif-quota",type=int,default=10)
    args=parser.parse_args()
    run(args.library,args.work,args.assets,args.step1,args.output,args.top_k,args.motif_quota)


if __name__=="__main__":
    main()
