"""Equal-size reference bootstrap gates; set scores are never per-peptide rewards."""
import numpy as np
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from .calibration import REGISTRY, calibrate, frozen_subsets, subseed
from .metric_scores import embedding_set_scores, sequence_set_scores

def evaluate_set(sequences, embeddings, reference_sequences, reference_embeddings,
                 anchor_embeddings, seed=42, max_n=2048, repeats=128, threshold=.5):
    if len(sequences) != len(embeddings) or len(reference_sequences) != len(reference_embeddings):
        raise ValueError("Sequence/embedding alignment mismatch")
    n = min(max_n, len(sequences), len(reference_sequences), len(anchor_embeddings))
    if n < 2:
        raise ValueError("Insufficient matched set size")
    anchors = np.asarray(anchor_embeddings, float)
    if not np.isfinite(anchors).all():
        raise ValueError("Nonfinite anchor embeddings")
    components = min(64, anchors.shape[1], len(anchors) - 1)
    pca = PCA(n_components=components, svd_solver="full")
    a = pca.fit_transform(anchors)
    ref = pca.transform(np.asarray(reference_embeddings, float))
    candidate = pca.transform(np.asarray(embeddings, float))
    anchor_ids = np.random.default_rng(subseed(seed, "set-anchors", n)).choice(len(a), n, replace=False)
    anchor_sample = a[anchor_ids]
    distances = pdist(anchor_sample, "sqeuclidean")
    bandwidth = float(np.median(distances[distances > 0]))
    if not np.isfinite(bandwidth):
        raise ValueError("Degenerate anchor embedding distribution")
    def scores(seqs, vectors):
        return {**sequence_set_scores(seqs),
                **embedding_set_scores(vectors, anchor_sample, bandwidth)}
    indices = frozen_subsets(len(reference_sequences), n, seed, repeats)
    reference_scores = {}
    for draw in indices:
        values = scores([reference_sequences[i] for i in draw], ref[draw])
        for key, value in values.items():
            reference_scores.setdefault(key, []).append(value)
    chosen = np.random.default_rng(subseed(seed, "set-candidate", n)).choice(
        len(sequences), n, replace=False)
    values = scores([sequences[i] for i in chosen], candidate[chosen])
    gate = calibrate(reference_scores, {key: [value] for key, value in values.items()},
                     1, threshold=threshold, level="set")
    return {"passed": bool(gate["passed"][0]), "sample_n": n, "repeats": repeats,
            "pca_components": components, "bandwidth_squared": bandwidth,
            "anchor_indices": anchor_ids.tolist(), "reference_indices": indices,
            "candidate_indices": chosen.tolist(), "reference_scores": reference_scores,
            "raw_scores": values,
            "percentiles": {name: float(p[0]) for name, p in gate["percentiles"].items()},
            "full_count": len(sequences), "full_unique_count": len(set(sequences))}
