"""Shared pure-sequence metric implementations for reference and generated peptides."""
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from rapidfuzz import process
from rapidfuzz.distance import Indel, Levenshtein
from .io import AA

def descriptors(sequences):
    if not sequences:
        raise ValueError("Empty candidate batch")
    for sequence in sequences:
        if not 8 <= len(sequence) <= 50:
            raise ValueError("Descriptors require length 8-50")
        if not set(sequence) <= set(AA):
            raise ValueError("Candidate contains non-canonical amino acids")
    rows = []
    for sequence in sequences:
        p = ProteinAnalysis(sequence)
        rows.append((len(sequence), p.charge_at_pH(7), p.gravy(), p.aromaticity(),
                     p.isoelectric_point(), p.molecular_weight()))
    names = ("length", "charge_ph7", "gravy", "aromaticity", "isoelectric_point", "molecular_weight")
    return {name: np.asarray(rows, dtype=float)[:, i] for i, name in enumerate(names)}

def anchor_novelty(sequences, anchors):
    if not anchors:
        raise ValueError("Empty anchor set")
    return np.asarray([1 - process.extractOne(
        sequence, anchors, scorer=Indel.normalized_similarity)[1] for sequence in sequences])

def cosine_nearest(queries, anchors, batch_size=128):
    queries, anchors = np.asarray(queries, float), np.asarray(anchors, float)
    if queries.ndim != 2 or anchors.ndim != 2 or queries.shape[1] != anchors.shape[1]:
        raise ValueError("Invalid embedding dimensions")
    if not np.isfinite(queries).all() or not np.isfinite(anchors).all():
        raise ValueError("Nonfinite embeddings")
    qnorm, anorm = np.linalg.norm(queries, axis=1), np.linalg.norm(anchors, axis=1)
    if not len(anchors) or np.any(qnorm == 0) or np.any(anorm == 0):
        raise ValueError("Empty or zero-norm embeddings")
    a = anchors / anorm[:, None]
    return np.concatenate([1 - np.clip((queries[i:i+batch_size] / qnorm[i:i+batch_size, None])
                                       @ a.T, -1, 1).max(axis=1)
                           for i in range(0, len(queries), batch_size)])

def sequence_set_scores(sequences):
    if len(sequences) < 2:
        raise ValueError("Set metrics require at least two sequences")
    distance = process.cdist(sequences, sequences, scorer=Levenshtein.normalized_distance,
                             workers=1, dtype=np.float64)
    return {"sequence_diversity": float(distance[np.triu_indices(len(sequences), 1)].mean()),
            "uniqueness": len(set(sequences)) / len(sequences)}

def embedding_set_scores(queries, anchors, bandwidth_squared):
    # Projection (if any) is fitted on anchors before calling this function.
    from scipy.spatial.distance import cdist
    q, a = np.asarray(queries, float), np.asarray(anchors, float)
    if min(len(q), len(a)) < 2 or not np.isfinite(q).all() or not np.isfinite(a).all():
        raise ValueError("Invalid set embeddings")
    if not np.isfinite(bandwidth_squared) or bandwidth_squared <= 0:
        raise ValueError("Bandwidth must be fixed from anchors and positive")
    cq, ca = np.atleast_2d(np.cov(q, rowvar=False)), np.atleast_2d(np.cov(a, rowvar=False))
    eigenvalues, vectors = np.linalg.eigh(np.atleast_2d(cq))
    root = (vectors * np.sqrt(np.maximum(eigenvalues, 0))) @ vectors.T
    cross = root @ np.atleast_2d(ca) @ root
    covariance_term = np.trace(cq) + np.trace(ca) - 2 * np.sqrt(
        np.maximum(np.linalg.eigvalsh((cross + cross.T) / 2), 0)).sum()
    fbd = float(np.sum((q.mean(0) - a.mean(0)) ** 2) + covariance_term)
    kernel = lambda x, y: np.exp(-cdist(x, y, "sqeuclidean") / (2 * bandwidth_squared)).mean()
    mmd = float(kernel(q, q) + kernel(a, a) - 2 * kernel(q, a))
    return {"embedding_fbd": max(0., fbd), "embedding_mmd": max(0., mmd)}
