"""Targeted contract tests: cold inference, unchanged gates and deterministic selection."""
import argparse
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from amp_submission.calibration import fingerprint
from amp_submission.generate import parser, too_similar
from amp_submission.hard_select import METRICS
from amp_submission.model_generate import merge_pool, sample_fresh, load_reference, choose
from amp_submission.io import digest


class FreshGenerationTests(unittest.TestCase):
    def test_default_is_full_mixed_entry(self):
        a = parser().parse_args([])
        self.assertEqual((a.arch, a.seed, a.n_sequences, a.top_k, a.pool_size, a.motif_quota),
                         ("mixed", 42, 50000, 100, 600000, 10))

    def test_console_entry_returns_success_not_report_as_exit_code(self):
        from amp_submission.generate import main
        with patch("sys.argv", ["generate"]), patch("amp_submission.model_generate.run", return_value={"ok": True}) as run:
            self.assertIsNone(main())
            run.assert_called_once()

    def test_merge_training_and_cross_branch_duplicates(self):
        known = "ACDEFGHIK"
        rows = merge_pool({"vq": [known, "KWVFKKLFK", "AAAAAAAAA"],
                           "dima": ["KWVFKKLFK"]}, {fingerprint(known)})
        self.assertFalse(rows[0]["hard_precheck"])
        self.assertEqual(rows[1]["source"], "dima|vq")
        self.assertEqual(rows[1]["raw_occurrences"], 2)
        self.assertFalse(rows[2]["hard_precheck"])

    def test_similarity_exact_boundary(self):
        self.assertFalse(too_similar("AAAAACCCCC", ["AAAAAGGCCC"]))
        self.assertTrue(too_similar("ACDEFGHIKL", ["ACDEFGHIKL"]))

    def test_manifest_tamper_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root / "reference-scores.json"
            path.write_text("{}")
            (root / "manifest.json").write_text(json.dumps({"sha256": {path.name: digest(path)}}))
            path.write_text('{"changed": 1}')
            with self.assertRaisesRegex(ValueError, "changed asset"):
                load_reference(root)

    def test_both_invocations_call_both_models_from_scratch(self):
        # Synthetic sampler only; this is not the full real-model certificate.
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            pdf = [0.0] * 41
            pdf[10] = 1.
            (root / "length-pdf.json").write_text(json.dumps(pdf))
            args = argparse.Namespace(assets=root, seed=42, threads=1, pool_size=256)
            def fake_load(assets, arch, device):
                from amp_submission.sampling import OBJECTIVES
                return {o: arch for o in OBJECTIVES}, None, {}
            def fake_sample(arch, policy, reference, bundle, lengths, device):
                return ["".join(np.random.choice(list("ACDEFGHIKLMNPQRSTVWY"), n)) for n in lengths]
            with patch("amp_submission.sampling.load_models", side_effect=fake_load) as loader, \
                 patch("amp_submission.sampling.sample", side_effect=fake_sample) as sampler, \
                 patch("torch.cuda.is_available", return_value=True):
                first = sample_fresh(args, root / "run1", lambda *a: None)
                # Existing unrelated output is deliberately not an input.
                second = sample_fresh(args, root / "run2", lambda *a: None)
            self.assertEqual(first, second)
            self.assertEqual(loader.call_count, 4)
            self.assertEqual(sampler.call_count, 4)
            self.assertEqual([c.args[1] for c in loader.call_args_list], ["vq", "dima", "vq", "dima"])
            for name in ("run1", "run2"):
                report = json.loads((root / name / "raw/complete.json").read_text())
                self.assertFalse(report["initial_library_used"])
                self.assertTrue(report["fresh_model_sampling"])

    def test_insufficient_pool_does_not_relax_or_fill(self):
        a = argparse.Namespace(n_sequences=50000)
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, "Insufficient qualified"):
                choose(pd.DataFrame({"sequence": ["ACDEFGHIK"]}),
                       np.ones((1, len(METRICS))), np.ones(1, bool), [], a, Path(temp))


if __name__ == "__main__":
    unittest.main()
