"""Length is a hard validity bound, never a reference-percentile quality target."""
import unittest
import numpy as np
from amp_submission.hard_select import METRICS, quality_arrays
from amp_submission.model_generate import merge_pool


class LengthHardOnlyTests(unittest.TestCase):
    def test_length_and_mass_are_not_quality_gates(self):
        self.assertNotIn("length", METRICS)
        self.assertNotIn("molecular_weight", METRICS)

    def test_length_and_mass_cannot_change_rank_or_pass(self):
        reference = {m: np.arange(1., 6.) for m in METRICS}
        reference.update(length=np.array([8., 9., 10., 20., 20.]),
                         molecular_weight=np.array([800., 900., 1000., 2000., 2000.]))
        candidates = {m: np.array([3., 3., 3.]) for m in METRICS}
        candidates.update(length=np.array([8., 40., 50.]),
                          molecular_weight=np.array([800., 4000., 5000.]))
        values, ranking, weakest = quality_arrays(reference, candidates)
        np.testing.assert_array_equal(values[0], values[1])
        np.testing.assert_array_equal(values[0], values[2])
        np.testing.assert_array_equal(ranking, np.repeat(ranking[0], 3))
        np.testing.assert_array_equal(weakest, np.repeat(weakest[0], 3))

    def test_generation_hard_length_accepts_41_to_50(self):
        sequences = [("ACDEFGHIKLMNPQRSTVWY" * 3)[:n] for n in (7, 8, 40, 41, 50, 51)]
        rows = merge_pool({"vq": sequences, "dima": []}, set())
        self.assertEqual([r["hard_precheck"] for r in rows],
                         [False, True, True, True, True, False])
