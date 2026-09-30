import unittest
import numpy as np
from amp_submission.hard_select import METRICS, quality_arrays, sequence_hard_mask, threshold_curve

class HardSelectTests(unittest.TestCase):
    def test_hard_sequence_boundaries(self):
        seqs = ["A"*7, "A"*8, "A"*50, "A"*51, "AXXXXXXX", "CCCCCCCC"]
        self.assertEqual(sequence_hard_mask(seqs, {"CCCCCCCC"}).tolist(),
                         [False, True, True, False, False, False])

    def test_threshold_and_missing(self):
        x = np.array([[.5, .5], [.3, .4], [.25, .2], [np.nan, 1.]])
        rows = threshold_curve(x, np.array([True, True, True, True]))
        self.assertEqual([r["eligible"] for r in rows], [1, 1, 1, 1, 2, 2])

    def test_hard_failure_cannot_be_rescued(self):
        rows = threshold_curve(np.ones((2, len(METRICS))), np.array([True, False]))
        self.assertTrue(all(r["eligible"] == 1 for r in rows))

    def test_partial_ranking_category_equal_weight(self):
        reference = {m: np.arange(1., 6.) for m in METRICS}
        candidates = {m: np.array([2., 3.]) for m in METRICS}
        values, ranking, weakest = quality_arrays(reference, candidates)
        np.testing.assert_allclose(ranking, (values[:, :-1].mean(axis=1) + values[:, -1]) / 2)
        np.testing.assert_allclose(weakest, values.min(axis=1))
        self.assertEqual(values.shape, (2, 5))

if __name__ == "__main__":
    unittest.main()
