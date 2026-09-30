"""Verified frozen model candidates and scores used by the default filtering entry."""
import json
import os
import tarfile
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from .calibration import fingerprint
from .hard_select import METRICS, THRESHOLDS, quality_arrays
from .io import digest, read_fasta, verify_assets
from .model_generate import merge_pool


def ensure_pool(root):
    root = Path(root)
    descriptor = Path("candidate-pool.json")
    release = json.loads(descriptor.read_text()) if descriptor.exists() else None
    if not root.exists():
        if release is None:
            raise RuntimeError("Frozen candidate pool is missing; provision the published pool")
        root.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".amp-download-", dir=root.parent) as temp:
            directory = Path(temp)
            archive = directory / "pool.tar.gz"
            print(json.dumps({"stage": "download/frozen_pool", "url": release["url"]}), flush=True)
            urllib.request.urlretrieve(release["url"], archive)
            if digest(archive) != release["sha256"]:
                raise ValueError("Frozen pool archive checksum mismatch")
            unpacked = directory / "unpacked"
            unpacked.mkdir()
            with tarfile.open(archive, "r:gz") as handle:
                for member in handle.getmembers():
                    target = (unpacked / member.name).resolve()
                    if (not member.isfile() or not target.is_relative_to(unpacked.resolve())):
                        raise ValueError("Unexpected frozen pool archive member")
                handle.extractall(unpacked, filter="data")
            verify_assets(unpacked)
            if digest(unpacked / "manifest.json") != release["manifest_sha256"]:
                raise ValueError("Frozen pool manifest checksum mismatch")
            os.replace(unpacked, root)
    manifest = verify_assets(root)
    if root == Path("frozen-pool") and release is not None:
        if digest(root / "manifest.json") != release["manifest_sha256"]:
            raise ValueError("Default pool differs from the published frozen pool")
    return manifest


def load_pool(root, manifest, reference):
    root = Path(root)
    if manifest.get("format") != "amp-frozen-pool-v1":
        raise ValueError("Unsupported frozen pool format")
    branches = {arch: read_fasta(root / (arch + ".fasta")) for arch in ("vq", "dima")}
    for arch, sequences in branches.items():
        if (len(sequences) != manifest["counts"][arch] or
                fingerprint(sequences) != manifest["raw_sha256"][arch]):
            raise ValueError("Frozen raw pool count/hash mismatch: " + arch)
    rows = merge_pool(branches, reference["training"])
    frame = pd.read_csv(root / "candidate-scores.csv", float_precision="round_trip")
    for field in ("sequence", "source", "raw_occurrences", "hard_precheck"):
        if frame[field].tolist() != [row[field] for row in rows]:
            raise ValueError("Frozen scores do not match raw candidates: " + field)
    percentiles, ranking, weakest = quality_arrays(reference["scores"], frame)
    possible = (frame.hard_precheck.to_numpy(bool) & np.isfinite(percentiles).all(axis=1)
                & (percentiles >= min(THRESHOLDS)).all(axis=1))
    audited = frame.known_amp_ratio_audited.to_numpy(bool)
    hard = frame.known_amp_ratio_pass.to_numpy(bool)
    if not np.array_equal(possible, audited) or np.any(hard & ~audited):
        raise ValueError("Frozen novelty audit coverage differs from the filtering protocol")
    frame["ranking"], frame["weakest"] = ranking, weakest
    for index, metric in enumerate(METRICS):
        frame[metric + "_percentile"] = percentiles[:, index]
    return branches, rows, frame, percentiles, hard
