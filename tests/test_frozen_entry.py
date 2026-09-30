"""Default entry must filter a frozen pool; model sampling requires an explicit flag."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from amp_submission.generate import main, parser


class FrozenEntryTests(unittest.TestCase):
    def test_default_does_not_request_sampling(self):
        args = parser().parse_args([])
        self.assertFalse(args.resample)
        self.assertEqual(args.candidate_pool, Path("frozen-pool"))

    def test_resampling_is_explicit(self):
        self.assertTrue(parser().parse_args(["--resample"]).resample)

    def test_default_missing_pool_fails_before_models(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("sys.argv", ["generate", "--candidate-pool", str(Path(temp) / "missing")]), \
                 patch("amp_submission.model_generate.sample_fresh") as sampler, \
                 patch("amp_submission.frozen_pool.urllib.request.urlretrieve", side_effect=FileNotFoundError):
                with self.assertRaises((FileNotFoundError, RuntimeError)):
                    main()
                sampler.assert_not_called()


if __name__ == "__main__":
    unittest.main()
