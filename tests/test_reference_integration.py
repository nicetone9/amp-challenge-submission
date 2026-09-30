import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from amp_submission.collection_gate import evaluate_set
from amp_submission.io import digest
from amp_submission.score_stage import load_scores

class ReferenceIntegrationTests(unittest.TestCase):
    def test_equal_sized_collection_reference_and_repeated_indices(self):
        seqs = ["ACDEFGHI", "KLMNPQRS", "YYYYYYYY", "AAAAAAAA", "KKLLKKLL", "GGLLGGLL"]
        embeddings = np.random.default_rng(42).normal(size=(6, 3))
        args = (seqs[:4], embeddings[:4], seqs, embeddings, embeddings + .1)
        a = evaluate_set(*args, repeats=4, max_n=3)
        self.assertEqual(a, evaluate_set(*args, repeats=4, max_n=3))
        self.assertEqual(a["sample_n"], 3)
        self.assertEqual(len(a["candidate_indices"]), 3)
        self.assertTrue(all(len(ids) == len(set(ids)) == 3 for ids in a["reference_indices"]))

    def test_manifest_length_column_matches_scored_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pd.DataFrame({"sequence": ["ACDEFGHI"], "length": [8]}).to_csv(
                root / "reference.csv", index=False)
            folder = root / "scores/descriptors"
            folder.mkdir(parents=True)
            for length in (8, 7):
                pd.DataFrame({"sequence": ["ACDEFGHI"], "length": [length]}).to_csv(
                    folder / "reference.csv", index=False)
                (folder / "complete.json").write_text(json.dumps(
                    {"sha256": {"reference": digest(folder / "reference.csv")}}))
                if length == 8:
                    result, _ = load_scores(root, "reference")
                    self.assertEqual(result.length.tolist(), [8])
                else:
                    with self.assertRaises(ValueError):
                        load_scores(root, "reference")

if __name__ == "__main__":
    unittest.main()
