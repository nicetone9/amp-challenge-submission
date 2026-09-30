"""Cluster-held-out exact 3-8mer enrichment against composition-preserving shuffles."""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from .calibration import fingerprint, subseed
from .io import digest
from .prepare_reference import save_json

def kmers(sequence):
    return {sequence[i:i+k] for k in range(3, 9) for i in range(len(sequence)-k+1)}

def bh(pvalues, total_tests=None):
    ordered = sorted(pvalues, key=lambda k: (pvalues[k], k))
    m = total_tests if total_tests is not None else len(ordered)
    result, running = {}, 1.
    for rank in range(len(ordered), 0, -1):
        name = ordered[rank-1]
        running = min(running, pvalues[name] * (m / rank))
        result[name] = running
    return result

def count_groups(frame, seed, label):
    positive, negative = Counter(), Counter()
    for sequence in frame.sequence:
        positive.update(kmers(sequence))
        rng = np.random.default_rng(subseed(seed, "motif-background", label, sequence))
        shuffled = "".join(rng.permutation(list(sequence)))
        negative.update(kmers(shuffled))
    return positive, negative

def enrichment(positive, negative, n, names, min_support, q_threshold, total_tests=None):
    pvalues, folds = {}, {}
    for motif in names:
        a, b = positive[motif], negative[motif]
        if a < min_support:
            continue
        folds[motif] = (a + .5) / (b + .5)
        pvalues[motif] = float(fisher_exact([[a, n-a], [b, n-b]], alternative="greater").pvalue)
    q = bh(pvalues, total_tests)
    return {motif: {"positive_support": positive[motif], "background_support": negative[motif],
                     "fold": folds[motif], "p": pvalues[motif], "q": q[motif]}
            for motif in sorted(pvalues) if q[motif] <= q_threshold and folds[motif] >= 2}

def discover(frame, seed=42):
    if not frame.sequence.is_unique or frame.cluster_id.isna().any():
        raise ValueError("Unique sequences and nonempty shared clusters required")
    clusters = sorted(frame.cluster_id.unique(), key=lambda c: (subseed(seed, "motif-split", c), c))
    if len(clusters) < 5:
        raise ValueError("Too few clusters for motif discovery/validation")
    discovery_clusters = set(clusters[:int(.8 * len(clusters))])
    discovery = frame[frame.cluster_id.isin(discovery_clusters)]
    validation = frame[~frame.cluster_id.isin(discovery_clusters)]
    pos, neg = count_groups(discovery, seed, "discovery")
    universe = set(pos) | set(neg)
    first = enrichment(pos, neg, len(discovery), universe, 10, .05, len(universe))
    vp, vn = count_groups(validation, seed, "validation")
    second = enrichment(vp, vn, len(validation), first, 5, .05, len(first))
    return {"seed": seed, "discovery_n": len(discovery), "validation_n": len(validation),
            "discovery_clusters": sorted(discovery_clusters),
            "validation_clusters": sorted(set(clusters) - discovery_clusters),
            "tested_discovery_universe": len(universe), "discovery_hits": len(first),
            "motifs": {m: {"discovery": first[m], "validation": second[m]} for m in second},
            "rules": {"widths": [3, 8], "q_max": .05, "fold_min": 2,
                      "discovery_min_support": 10, "validation_min_support": 5,
                      "background": "one deterministic composition-preserving shuffle per sequence"},
            "warning": "Independent supporting evidence only; not proof of AMP activity or a qualification gate."}

def scan(sequence, motifs):
    return [{"motif": sequence[i:i+k], "start": i, "end": i+k}
            for i in range(len(sequence)) for k in range(3, 9)
            if sequence[i:i+k] in motifs and i+k <= len(sequence)]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("work/six-metrics"))
    p.add_argument("--seed", type=int, default=42)
    a = p.parse_args()
    report = discover(pd.read_csv(a.work / "reference.csv"), a.seed)
    report["reference_sha256"] = digest(a.work / "reference.csv")
    report["code_sha256"] = digest(__file__)
    target = a.work / "motifs"
    save_json(target / "frozen.json", report)
    candidates = pd.read_csv(a.work / "candidates.csv", keep_default_na=False)
    motifs = report["motifs"]
    records = [{"sequence_sha256": fingerprint(s), "hits": scan(s, motifs)}
               for s in candidates.sequence]
    save_json(target / "candidate-hits.json", records)
    summary = {"status": "SUPPORT_ONLY", "reliable_motifs": len(motifs),
               "candidate_hit_count": sum(bool(row["hits"]) for row in records),
               "candidate_count": len(records), "frozen_sha256": digest(target / "frozen.json")}
    save_json(target / "summary.json", summary)
    print(summary)

if __name__ == "__main__":
    main()
