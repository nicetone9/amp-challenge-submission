"""Frozen empirical calibration; no experimental labels enter predicted-score CDFs."""
import hashlib
import json
from dataclasses import asdict, dataclass
import numpy as np

CATEGORIES = ("physicochemical", "potency", "embedding", "synthesizability", "novelty", "diversity")

@dataclass(frozen=True)
class Metric:
    name: str
    category: str
    direction: str
    level: str
    unit: str
    implementation: str
    status: str = "IMPLEMENTED"
    reason: str = ""

REGISTRY = tuple(
    Metric(name, "physicochemical", "typical", "sequence", unit, implementation)
    for name, unit, implementation in (
        ("length", "aa", "len; complete sequence"),
        ("charge_ph7", "e", "Bio.SeqUtils.ProtParam.ProteinAnalysis.charge_at_pH(7.0)"),
        ("gravy", "Kyte-Doolittle", "Bio.SeqUtils.ProtParam.ProteinAnalysis.gravy()"),
        ("aromaticity", "fraction", "Bio.SeqUtils.ProtParam.ProteinAnalysis.aromaticity()"),
        ("isoelectric_point", "pH", "Bio.SeqUtils.ProtParam.ProteinAnalysis.isoelectric_point()"),
        ("molecular_weight", "Da", "Bio.SeqUtils.ProtParam.ProteinAnalysis.molecular_weight(); free termini"),
    )
) + tuple(Metric("ania_" + species, "potency", "low", "sequence",
                 "ANIA native Predicted Log MIC Value", "ANIA weights and code bound by asset manifest")
          for species in ("ecoli", "paeruginosa", "saureus")) + (
    Metric("hemopi2_hc50_um", "potency", "high", "sequence", "uM",
           "HemoPI2 regression; safety adjunct, not MIC or measured HC50"),
    Metric("embedding_nn_cosine_distance", "embedding", "low", "sequence", "cosine distance",
           "frozen ProtT5 mean residue embedding; nearest disjoint AMP anchor"),
    Metric("peptide_synthesis_score", "synthesizability", "high", "sequence", "UNAVAILABLE",
           "UNAVAILABLE", "UNAVAILABLE", "No validated short-peptide synthesis implementation/weights selected"),
    Metric("anchor_novelty", "novelty", "high", "sequence", "1 - Levenshtein.ratio",
           "rapidfuzz Indel.normalized_similarity; disjoint AMP anchors"),
    Metric("embedding_fbd", "embedding", "low", "set", "squared embedding distance",
           "Gaussian Frechet distance; PCA64 fit only on anchors"),
    Metric("embedding_mmd", "embedding", "low", "set", "RBF MMD squared, biased",
           "fixed anchor median-distance bandwidth; equal sample size"),
    Metric("sequence_diversity", "diversity", "high", "set", "mean normalized edit distance",
           "rapidfuzz Levenshtein.normalized_distance; unordered distinct pairs"),
    Metric("uniqueness", "diversity", "high", "set", "fraction", "unique complete sequences / count"),
)

def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()

def subseed(seed, *parts):
    return int(fingerprint([int(seed), *parts])[:16], 16)

def registry_payload():
    return [asdict(metric) for metric in REGISTRY]

def midpoint_cdf(reference, values):
    ref = np.asarray(reference, dtype=float)
    if ref.ndim != 1 or not len(ref) or not np.isfinite(ref).all():
        raise ValueError("Reference must be a nonempty finite vector")
    ref = np.sort(ref)
    x = np.asarray(values, dtype=float)
    cdf = (np.searchsorted(ref, x, side="left") +
           np.searchsorted(ref, x, side="right")) / (2 * len(ref))
    return np.where(np.isfinite(x), cdf, np.nan)

def percentile(reference, values, direction):
    p = midpoint_cdf(reference, values)
    if direction == "high":
        return p
    if direction == "low":
        return 1 - p
    if direction == "typical":
        reference_p = midpoint_cdf(reference, reference)
        typical_reference = 2 * np.minimum(reference_p, 1 - reference_p)
        typical_values = 2 * np.minimum(p, 1 - p)
        return midpoint_cdf(typical_reference, typical_values)
    raise ValueError("Unknown direction: " + direction)

def split_clusters(cluster_ids, seed):
    clusters = sorted(set(cluster_ids), key=lambda c: (subseed(seed, "anchor", c), c))
    if len(clusters) < 2 or any(not c for c in clusters):
        raise ValueError("At least two nonempty cluster IDs required")
    anchor = set(clusters[:len(clusters) // 2])
    return ["anchor" if c in anchor else "calibration" for c in cluster_ids]

def frozen_subsets(population_size, requested_size, seed, repeats=128):
    n = min(population_size, requested_size)
    if n < 2 or repeats < 1:
        raise ValueError("Insufficient reference for a set-level comparison")
    return [np.random.default_rng(subseed(seed, "set-reference", n, i)).choice(
        population_size, n, replace=False).tolist() for i in range(repeats)]

def calibrate(reference_scores, candidate_scores, size, threshold=.5, level="sequence", registry=REGISTRY):
    if not 0 <= threshold <= 1:
        raise ValueError("Invalid percentile threshold")
    percentiles, statuses, categories = {}, {}, {}
    for metric in registry:
        if metric.level != level:
            continue
        name = metric.name
        reference = reference_scores.get(name)
        values = candidate_scores.get(name)
        if metric.status != "IMPLEMENTED" or reference is None or values is None:
            p = np.full(size, np.nan)
            statuses[name] = metric.reason or "Missing reference or candidate scores"
        else:
            values = np.asarray(values, dtype=float)
            if values.shape != (size,):
                raise ValueError("Score length mismatch: " + name)
            try:
                p = percentile(reference, values, metric.direction)
                statuses[name] = "COMPUTED"
            except ValueError as error:
                p = np.full(size, np.nan)
                statuses[name] = str(error)
        percentiles[name] = p
        categories.setdefault(metric.category, []).append(p)
    if not percentiles:
        raise ValueError("No registered metrics for requested level")
    matrix = np.asarray(list(percentiles.values()))
    passed = (np.isfinite(matrix) & (matrix >= threshold)).all(axis=0)
    category_scores = {key: np.mean(values, axis=0) for key, values in categories.items()}
    # Missing metrics propagate NaN; never rank by a silently reduced objective.
    ranking = np.mean(list(category_scores.values()), axis=0)
    return dict(percentiles=percentiles, statuses=statuses, passed=passed,
                category_scores=category_scores, ranking=ranking, weakest=np.min(matrix, axis=0))
