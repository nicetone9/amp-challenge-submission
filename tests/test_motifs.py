import unittest
import pandas as pd
from amp_submission.motifs import bh, discover, kmers, scan

class MotifTests(unittest.TestCase):
    def test_bh_correction(self):
        self.assertEqual(bh({"a": .01, "b": .04, "c": .05}),
                         {"c": .05, "b": .05, "a": .03})

    def test_exact_positions(self):
        self.assertEqual(scan("ACDACDAC", {"ACD"}),
                         [{"motif": "ACD", "start": 0, "end": 3},
                          {"motif": "ACD", "start": 3, "end": 6}])
        self.assertTrue(all(3 <= len(x) <= 8 for x in kmers("ACDEFGHIKL")))

    def test_cluster_discovery_validation_disjoint(self):
        frame = pd.DataFrame({"sequence": ["ACDEFGHI", "CDEFGHIK", "DEFGHIKL",
                                          "EFGHIKLM", "FGHIKLMN", "GHIKLMNP"],
                              "cluster_id": ["a", "a", "b", "c", "d", "e"]})
        a = discover(frame)
        self.assertEqual(a, discover(frame.iloc[::-1]))
        self.assertFalse(set(a["discovery_clusters"]) & set(a["validation_clusters"]))
        self.assertEqual(a["discovery_n"] + a["validation_n"], len(frame))
        self.assertEqual(a["motifs"], {})

if __name__ == "__main__":
    unittest.main()
