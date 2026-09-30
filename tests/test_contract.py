import json
import tempfile
import unittest
from pathlib import Path
import Levenshtein
from amp_submission.generate import parser, too_similar, rank_top, validate
from amp_submission.io import digest, verify_assets, write_fasta, read_fasta
from amp_submission.scoring import rewards

class ContractTests(unittest.TestCase):
    def test_defaults(self):
        a = parser().parse_args([])
        self.assertEqual((a.arch, a.seed, str(a.output)), ("mixed", 42, "generate"))

    def test_exact_boundary(self):
        for a, b in [("AAAA", "AAA"), ("AAAAA", "AAA"), ("ACDEFGHIKL", "ACDEFGHILL")]:
            self.assertEqual(too_similar(a, [b]), Levenshtein.ratio(a, b) > .8)

    def test_stable_ranking(self):
        pool = ["ACDEFGHI", "KKKKRRRR", "WWWWYYYY"]
        self.assertEqual(rank_top(pool, [.5, .5, .5], 2), rank_top(pool, [.5, .5, .5], 2))
        self.assertEqual(rank_top(pool, [.5, .5, .5], 1), pool[:1])

    def test_repeated_fasta_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "library.fasta"
            seqs = ["ACDEFGHI", "KLMNPQRS"]
            write_fasta(p, seqs)
            before = p.read_bytes()
            write_fasta(p, seqs)
            self.assertEqual(before, p.read_bytes())
            self.assertEqual(read_fasta(p), seqs)

    def test_validation(self):
        seqs = ["ACDEFGHI", "KLMNPQRS"]
        validate(seqs, seqs[:1], ["YYYYYYYY"], 2, 1)
        for bad in [seqs[:1] * 2, ["SHORT", "KLMNPQRS"], ["ACDEFGHX", "KLMNPQRS"]]:
            with self.assertRaises(ValueError):
                validate(bad, bad[:1], ["YYYYYYYY"], 2, 1)
        with self.assertRaises(ValueError):
            validate(seqs, seqs[:1], seqs[:1], 2, 1)

    def test_asset_integrity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "x").write_text("good")
            (root / "manifest.json").write_text(json.dumps({"sha256": {"x": digest(root / "x")}}))
            verify_assets(root)
            (root / "x").write_text("bad")
            with self.assertRaises(ValueError):
                verify_assets(root)

    def test_reward_direction_and_joint(self):
        import numpy as np
        keys = ["ania_ecoli", "ania_paeruginosa", "ania_saureus", "hemopi2_hc50_um"]
        cal = {k: [1, 2, 3] for k in keys}
        scores = {k: np.array([1., 3.]) for k in keys}
        result = rewards(scores, cal)
        self.assertGreater(result["broad"][0], result["broad"][1])
        self.assertTrue(np.allclose(result["joint"],
            np.mean([result[k] for k in ("broad", "gram_positive", "gram_negative", "mdr", "selectivity")], axis=0)))

if __name__ == "__main__":
    unittest.main()
