import tempfile
import unittest
from pathlib import Path

import numpy as np

from amp_submission.hard_audit import audit_batch
from amp_submission.score_stage import cache_scores


class HardAuditTests(unittest.TestCase):
    def test_boundary_exact_hit_and_novel_sequence(self):
        # Equal-length one substitution gives ratio exactly 0.8.
        result = audit_batch(["ACDEF", "ACDEG", "WWWWW"], ["ACDEF"])
        np.testing.assert_array_equal(result["official_reference_pass"], [0, 1, 1])

    def test_resume_preserves_order_and_results(self):
        with tempfile.TemporaryDirectory() as directory:
            calls = []
            def compute(sequences):
                calls.append(tuple(sequences))
                return audit_batch(sequences, ["ACDEF"])
            sequences = ["WWWWW", "ACDEF", "ACDEG"]
            cold = cache_scores(sequences, Path(directory), {"reference": "fixed"}, compute, 2)
            warm = cache_scores(sequences, Path(directory), {"reference": "fixed"}, compute, 2)
            self.assertEqual(len(calls), 2)
            np.testing.assert_array_equal(cold["official_reference_pass"], warm["official_reference_pass"])


if __name__ == "__main__":
    unittest.main()
