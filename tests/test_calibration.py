import unittest
import numpy as np
from amp_submission.calibration import (
    CATEGORIES, REGISTRY, Metric, calibrate, fingerprint, frozen_subsets,
    midpoint_cdf, percentile, split_clusters, subseed,
)
from amp_submission.selection import allocate, assign_tiers, select

class CalibrationTests(unittest.TestCase):
    def test_six_categories_required(self):
        self.assertEqual(set(CATEGORIES), {m.category for m in REGISTRY})
        self.assertTrue(any(m.category == "novelty" and m.level == "sequence" for m in REGISTRY))

    def test_midpoint_and_fifty_boundary(self):
        np.testing.assert_allclose(midpoint_cdf([1, 2, 2, 3], [0, 1, 2, 3, 4]),
                                   [0, .125, .5, .875, 1])
        metric = Metric("x", "potency", "high", "sequence", "x", "test")
        out = calibrate({"x": [1, 2, 2, 3]}, {"x": [2]}, 1, registry=[metric])
        self.assertTrue(out["passed"][0])

    def test_directions(self):
        np.testing.assert_allclose(percentile([1, 2, 3], [1, 3], "low"), [5/6, 1/6])
        np.testing.assert_allclose(percentile([1, 2, 3], [1, 3], "high"), [1/6, 5/6])

    def test_typicality_not_monotone(self):
        p = percentile([0, 1, 2, 3, 4], [-1, 2, 5], "typical")
        self.assertGreater(p[1], p[0])
        self.assertGreater(p[1], p[2])
        self.assertTrue(np.isnan(percentile([1, 2], [np.nan], "typical")[0]))

    def test_nonfinite_rejected(self):
        for ref in ([], [1, np.nan], [np.inf]):
            with self.assertRaises(ValueError):
                midpoint_cdf(ref, [1])
        for x in (np.nan, np.inf, -np.inf):
            self.assertTrue(np.isnan(midpoint_cdf([1, 2], [x])[0]))

    def test_and_not_average(self):
        metrics = [Metric(x, "potency", "high", "sequence", "", "") for x in ("x", "y")]
        out = calibrate({"x": [0, 1], "y": [0, 1]}, {"x": [2], "y": [-1]}, 1, registry=metrics)
        self.assertEqual(out["ranking"][0], .5)
        self.assertFalse(out["passed"][0])

    def test_missing_blocks(self):
        out = calibrate({}, {}, 3)
        self.assertFalse(out["passed"].any())
        self.assertTrue(np.isnan(out["ranking"]).all())
        self.assertIn("peptide_synthesis_score", out["statuses"])

    def test_equal_category_weights(self):
        metrics = [Metric(x, c, "high", "sequence", "", "") for x, c in
                   (("a", "potency"), ("b", "potency"), ("c", "novelty"))]
        out = calibrate({k: [0, 1] for k in "abc"}, {"a": [2], "b": [2], "c": [-1]}, 1, registry=metrics)
        self.assertEqual(out["ranking"][0], .5)

    def test_cluster_split_disjoint_and_order_invariant(self):
        ids = ["a", "a", "b", "c", "d"]
        parts = split_clusters(ids, 42)
        self.assertEqual(parts[0], parts[1])
        self.assertEqual(parts, list(reversed(split_clusters(list(reversed(ids)), 42))))
        self.assertFalse({c for c, p in zip(ids, parts) if p == "anchor"} &
                         {c for c, p in zip(ids, parts) if p == "calibration"})

    def test_frozen_subsets(self):
        draws = frozen_subsets(23, 100, 42)
        self.assertEqual(len(draws), 128)
        self.assertEqual(draws, frozen_subsets(23, 100, 42))
        self.assertTrue(all(len(x) == len(set(x)) == 23 for x in draws))
        with self.assertRaises(ValueError):
            frozen_subsets(1, 100, 42)

    def test_namespace_hashes(self):
        self.assertEqual(fingerprint({"a": 1, "b": 2}), fingerprint({"b": 2, "a": 1}))
        self.assertNotEqual(subseed(42, "vq"), subseed(42, "dima"))

class SelectionTests(unittest.TestCase):
    def test_capacity_and_rounding(self):
        self.assertEqual(allocate({"a": 1, "b": 9}, 9), {"a": 1, "b": 8})
        self.assertEqual(allocate({"b": 1, "a": 1}, 1), {"a": 1, "b": 0})
        self.assertEqual(sum(allocate({"a": 3, "b": 17, "c": 1}, 13).values()), 13)
        with self.assertRaises(ValueError):
            allocate({"a": 1}, 2)

    def rows(self):
        return [dict(sequence=f"seq{i}", ranking=.5 + i/100, weakest=.5, passed=True,
                     cluster=str(i//6), motif_support=0) for i in range(12)]

    def test_high_low_odd_and_top_membership(self):
        rows = assign_tiers(self.rows(), 42)
        library = select(rows, 9)
        self.assertEqual(sum(r["tier"] == "high" for r in library), 5)
        top = select(library, 5)
        self.assertEqual(sum(r["tier"] == "high" for r in top), 3)
        self.assertTrue({r["sequence"] for r in top} <= {r["sequence"] for r in library})
        self.assertEqual(library, select(list(reversed(rows)), 9))

    def test_low_still_qualified(self):
        rows = self.rows()
        rows[0]["passed"] = False
        with self.assertRaises(ValueError):
            assign_tiers(rows, 42)

    def test_singleton_determinism(self):
        rows = [dict(r, cluster=r["sequence"]) for r in self.rows()]
        self.assertEqual(assign_tiers(rows, 42), assign_tiers(list(reversed(rows)), 42))

    def test_no_quota_relaxation(self):
        rows = [dict(r, tier="high") for r in self.rows()]
        with self.assertRaises(ValueError):
            select(rows, 4)

if __name__ == "__main__":
    unittest.main()
