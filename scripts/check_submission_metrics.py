"""Create a submission-file checklist and descriptive metrics table from verified output."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd

from amp_submission.calibration import fingerprint
from amp_submission.hard_select import METRICS, quality_arrays
from amp_submission.io import AA, digest, read_fasta
from amp_submission.model_generate import load_reference
from amp_submission.prepare_reference import save_json


def check(delivery, reference_root, output, verification_path=None):
    delivery, output = Path(delivery), Path(output)
    report_path = delivery / "generation.json"
    if not report_path.exists():
        report_path = delivery / "complete.json"
    report = json.loads(report_path.read_text())
    certificate = json.loads(Path(verification_path).read_text()) if verification_path else None
    reference = load_reference(reference_root)
    legacy = "reference_sha256" not in report
    if legacy:
        if (not certificate or certificate["fasta_sha256"] != report["fasta_sha256"]
                or not certificate["byte_identical"] or not certificate["whole_library_known_reference_ratio_le_0_8"]
                or not certificate["independent_percentile_and_rank_recalculation"]):
            raise ValueError("Legacy delivery requires its verified unchanged-FASTA certificate")
    elif report["reference_sha256"] != reference["manifest_sha256"]:
        raise ValueError("Selection reference differs from generation")
    frames, sequences = {}, {}
    for name, expected in (("library", 50000), ("top", 100)):
        fasta = delivery / (name + ".fasta")
        if digest(fasta) != report["fasta_sha256"][name]:
            raise ValueError("Delivery FASTA changed: " + name)
        sequences[name] = read_fasta(fasta)
        frame = pd.read_csv(delivery / (name + "-scores.csv"), float_precision="round_trip")
        if frame.sequence.tolist() != sequences[name]:
            raise ValueError("Score table differs from exported FASTA")
        if len(sequences[name]) != expected:
            raise ValueError("Wrong submission quantity")
        frames[name] = frame
    threshold = report.get("threshold", report.get("reference_percentile"))
    rows = []
    def add(metric, standard, counts, ranges=None, state="PASS"):
        rows.append({"metric": metric, "requirement": standard,
                     "library": str(counts[0]), "top100": str(counts[1]),
                     "library_range": ranges[0] if ranges else "",
                     "top100_range": ranges[1] if ranges else "", "status": state})
    sizes = [len(sequences[name]) for name in ("library", "top")]
    add("数量", "library=50000；Top100=100", sizes)
    add("唯一序列", "无重复", [len(set(sequences[name])) for name in ("library", "top")])
    add("标准氨基酸", "仅20种标准AA", [sum(set(s) <= set(AA) for s in sequences[name]) for name in ("library", "top")])
    lengths = [pd.Series([len(s) for s in sequences[name]]) for name in ("library", "top")]
    add("长度合法性", "8–50 aa；仅硬性合法范围", [int(v.between(8, 50).sum()) for v in lengths],
        [f"{v.min()}–{v.max()} aa" for v in lengths])
    add("低复杂度排除", "单一AA占比≤0.6",
        [sum(max(Counter(s).values()) / len(s) <= .6 for s in sequences[name]) for name in ("library", "top")])
    add("训练集精确重叠排除", "训练序列哈希不命中",
        [sum(fingerprint(s) not in reference["training"] for s in sequences[name]) for name in ("library", "top")])
    add("已知AMP相似性", "Levenshtein ratio≤0.8；生成时完整库审核已通过",
        [int(frames[name].known_amp_ratio_pass.sum()) if not legacy else len(sequences[name]) for name in ("library", "top")])
    if not legacy:
        command = report["cd_hit_command"]
        for flag, value in (("-c", "0.5"), ("-aS", "0.8"), ("-aL", "0.8"), ("-T", "1")):
            if command[command.index(flag) + 1] != value:
                raise ValueError("Clustering protocol changed")
    add("CD-HIT聚类", "identity=0.50；双向coverage=0.80；单线程",
        [f"{frames[name].cluster.nunique()} clusters" for name in ("library", "top")])
    add("Top100属于library", "全部来自本次library", [sizes[0], len(set(sequences["top"]) & set(sequences["library"]))])
    tiers = [Counter(frames[name].tier) for name in ("library", "top")]
    add("高/低组配额", "library各25000；Top100各50",
        [f"high={value['high']}; low={value['low']}" for value in tiers])
    add("短motif支持", "Top100至少10；3–4 aa motif；每组至少5个配额位",
        [int(frames[name].short_motif_support.gt(0).sum()) for name in ("library", "top")])
    percentile_values = {}
    for name in ("library", "top"):
        values = (frames[name][[metric + "_percentile" for metric in METRICS]].to_numpy(float) if legacy
                  else quality_arrays(reference["scores"], frames[name])[0])
        percentile_values[name] = values
        if not np.isfinite(values).all() or not (values >= threshold).all():
            raise ValueError("Exported sequences fail the frozen quality threshold")
    for index, metric in enumerate(METRICS):
        minimums = [percentile_values[name][:, index].min() for name in ("library", "top")]
        ranges = []
        for name in ("library", "top"):
            values = frames[name][metric + "_percentile"] if legacy else frames[name][metric]
            label = "percentile " if legacy else ""
            ranges.append(f"{label}{values.min():.6g}–{values.max():.6g}; median={values.median():.6g}")
        add(metric, f"校准percentile≥{threshold:.2f}；最小值 library={minimums[0]:.6g}, Top={minimums[1]:.6g}",
            sizes, ranges)
    add("分子量", "仅描述；不参与质量门槛或排名", ["描述项", "描述项"],
        [f"{frames[name].molecular_weight.min():.6g}–{frames[name].molecular_weight.max():.6g} Da" for name in ("library", "top")],
        state="DESCRIPTIVE")
    for metric in ("抗菌活性/MIC与溶血安全性", "ProtT5/集合表征指标", "合成可行性", "湿实验活性与安全性"):
        add(metric, "当前输出未完成此项验证", ["未验证", "未验证"], state="UNVERIFIED")
    motif_hits = frames["top"].short_motif_support.gt(0)
    if motif_hits.sum() < 10 or any((motif_hits & frames["top"].tier.eq(tier)).sum() < 5 for tier in ("high", "low")):
        raise ValueError("Short-motif quota failed")
    for row in rows[:9]:
        if (row["library"].isdigit() and int(row["library"]) != sizes[0]
                or row["top100"].isdigit() and int(row["top100"]) != sizes[1]):
            raise ValueError("Hard checklist failed: " + row["metric"])
    if (tiers[0] != Counter(high=25000, low=25000) or tiers[1] != Counter(high=50, low=50)
            or not set(sequences["top"]) <= set(sequences["library"])):
        raise ValueError("Tier or Top100 membership check failed")
    output.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    table.to_csv(output / "submission-metrics-check.csv", index=False)
    header = "| 指标 | 标准 | 5万库 | Top100 | 5万库范围 | Top100范围 | 状态 |\n|---|---|---:|---:|---|---|---|\n"
    markdown = header + "".join("| " + " | ".join(str(value).replace("|", "/") for value in row.values()) + " |\n" for row in rows)
    (output / "submission-metrics-check.md").write_text(markdown)
    summary = {"fasta_sha256": report["fasta_sha256"], "threshold": threshold,
               "quality_metrics": list(METRICS), "library_count": sizes[0], "top_count": sizes[1],
               "hard_checks_passed": True, "all_six_categories_passed": False,
               "reference_manifest_sha256": reference["manifest_sha256"],
               "delivery_report_sha256": digest(report_path),
               "reused_verification_sha256": digest(verification_path) if verification_path else None,
               "novelty_evidence": "The unchanged FASTAs passed generation.validate against the frozen combined known-AMP reference",
               "rows": rows}
    save_json(output / "submission-metrics-check.json", summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path("checkpoint/generation-reference"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verification", type=Path, help="Existing unchanged-FASTA certificate for a legacy fixed-pool delivery")
    args = parser.parse_args()
    check(args.delivery, args.reference, args.output, args.verification)
