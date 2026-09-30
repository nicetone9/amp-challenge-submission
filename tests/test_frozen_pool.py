"""Synthetic input-integrity tests; full library validation is performed separately."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from amp_submission.calibration import fingerprint
from amp_submission.frozen_pool import load_pool
from amp_submission.hard_select import METRICS
from amp_submission.io import write_fasta
from amp_submission.model_generate import merge_pool


class FrozenPoolIntegrityTests(unittest.TestCase):
    def fixture(self, root):
        branches = {"vq": ["ACDEFGHIK"], "dima": ["LMNPQRSTVW"]}
        for arch, sequences in branches.items():
            write_fasta(root / (arch + ".fasta"), sequences)
        reference = {"training": set(), "scores": {metric: [0., 1., 2.] for metric in METRICS}}
        rows = merge_pool(branches, reference["training"])
        frame = pd.DataFrame(rows)
        for metric in METRICS:
            frame[metric] = 1.
        frame["known_amp_ratio_audited"] = True
        frame["known_amp_ratio_pass"] = True
        frame.to_csv(root / "candidate-scores.csv", index=False)
        manifest = {"format": "amp-frozen-pool-v1", "counts": {"vq": 1, "dima": 1},
                    "raw_sha256": {arch: fingerprint(sequences) for arch, sequences in branches.items()}}
        return manifest, reference, frame

    def test_valid_pool_keeps_raw_order_and_recomputes_quality(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, reference, _ = self.fixture(root)
            branches, rows, frame, percentiles, hard = load_pool(root, manifest, reference)
            self.assertEqual(frame.sequence.tolist(), branches["vq"] + branches["dima"])
            self.assertEqual(len(rows), 2)
            self.assertTrue((percentiles >= .25).all())
            self.assertTrue(hard.all())

    def test_changed_raw_sequences_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, reference, _ = self.fixture(root)
            write_fasta(root / "vq.fasta", ["ACDEFGHIL"])
            with self.assertRaisesRegex(ValueError, "count/hash mismatch"):
                load_pool(root, manifest, reference)

    def test_scores_from_another_pool_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, reference, frame = self.fixture(root)
            frame["sequence"] = frame.sequence.tolist()[::-1]
            frame.to_csv(root / "candidate-scores.csv", index=False)
            with self.assertRaisesRegex(ValueError, "do not match raw candidates"):
                load_pool(root, manifest, reference)

    def test_selectable_but_unaudited_candidates_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, reference, frame = self.fixture(root)
            frame.loc[0, "known_amp_ratio_audited"] = False
            frame.to_csv(root / "candidate-scores.csv", index=False)
            with self.assertRaisesRegex(ValueError, "novelty audit coverage"):
                load_pool(root, manifest, reference)


if __name__ == "__main__":
    unittest.main()
