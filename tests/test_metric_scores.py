import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from amp_submission.metric_scores import descriptors, anchor_novelty, cosine_nearest, sequence_set_scores, embedding_set_scores
from amp_submission.score_stage import cache_scores
from amp_submission.selection import cluster_sequences

class MetricTests(unittest.TestCase):
    def test_shared_descriptor_values(self):
        seqs = ["ACDEFGHI", "KKLLKKLL"]
        scores = descriptors(seqs)
        np.testing.assert_equal(scores["length"], [8, 8])
        np.testing.assert_allclose(scores["aromaticity"], [1/8, 0])
        for key in scores:
            np.testing.assert_equal(scores[key], descriptors(seqs)[key])
        with self.assertRaises(ValueError):
            descriptors(["A" * 41])

    def test_novelty_direction_and_cosine(self):
        p = anchor_novelty(["ACDEFGHI", "YYYYYYYY"], ["ACDEFGHI"])
        self.assertEqual(p[0], 0)
        self.assertGreater(p[1], p[0])
        np.testing.assert_allclose(cosine_nearest([[1, 0], [0, 1]], [[1, 0]]), [0, 1])

    def test_set_metrics_symmetry_and_identity(self):
        a = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]])
        out = embedding_set_scores(a, a, 1.)
        self.assertLess(out["embedding_fbd"], 1e-10)
        self.assertLess(out["embedding_mmd"], 1e-10)
        b = a + 1
        self.assertEqual(embedding_set_scores(a, b, 1.), embedding_set_scores(b, a, 1.))
        self.assertEqual(sequence_set_scores(["AAAAAAAA", "AAAAAAAA"])["sequence_diversity"], 0)
        self.assertEqual(sequence_set_scores(["AAAAAAAA", "AAAAAAAA"])["uniqueness"], .5)

    def test_cache_warm_cold_resume_same_and_no_rng_effect(self):
        with tempfile.TemporaryDirectory() as tmp:
            seqs = ["ACDEFGHI", "KKLLKKLL", "YYYYYYYY"]
            np.random.seed(42)
            state = np.random.get_state()
            a = cache_scores(seqs, tmp, {"code": "a"}, descriptors, 2)
            self.assertEqual(np.random.get_state()[2:], state[2:])
            np.testing.assert_equal(np.random.get_state()[1], state[1])
            b = cache_scores(seqs, tmp, {"code": "a"},
                             lambda _: self.fail("warm cache recomputed"), 2)
            for key in a:
                np.testing.assert_equal(a[key], b[key])
            # Interrupted after the first batch: already stored batch then resume.
            with tempfile.TemporaryDirectory() as partial:
                cache_scores(seqs[:2], partial, {"code": "a"}, descriptors, 2)
                c = cache_scores(seqs, partial, {"code": "a"}, descriptors, 2)
                for key in a:
                    np.testing.assert_equal(a[key], c[key])

    def test_cache_corruption_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            cache_scores(["ACDEFGHI"], tmp, {}, descriptors)
            p = next(Path(tmp).glob("*.json"))
            value = json.loads(p.read_text())
            value["scores"]["length"] = [123]
            p.write_text(json.dumps(value))
            with self.assertRaises(ValueError):
                cache_scores(["ACDEFGHI"], tmp, {}, descriptors)

    def test_cdhit_keeps_short_peptides_and_is_order_invariant(self):
        seqs = ["ACDEFGHI", "ACDEFGHIK", "ACDEFGHIKL", "YYYYYYYY", "KKLLKKLL"]
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            first, cmd = cluster_sequences(seqs, a)
            second, _ = cluster_sequences(list(reversed(seqs)), b)
            self.assertEqual(set(first), set(seqs))
            self.assertEqual(first, second)
            self.assertEqual(cmd[cmd.index("-l")+1], "7")
            self.assertEqual(cmd[cmd.index("-T")+1], "1")

if __name__ == "__main__":
    unittest.main()
