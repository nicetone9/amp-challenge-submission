import tempfile
import unittest
from pathlib import Path
import numpy as np
from amp_submission.expand_pool import extend_batches

class PoolResumeTests(unittest.TestCase):
    def test_resume_cache_budget_extension_identical_prefix(self):
        def sample(arch, index, n, seed):
            rng = np.random.default_rng(seed)
            return ["".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 8)) for _ in range(n)]
        with tempfile.TemporaryDirectory() as cold, tempfile.TemporaryDirectory() as resumed:
            identity = {"seed": 42, "model": "fixture"}
            extend_batches(cold, identity, 10, sample, 4)
            extend_batches(resumed, identity, 3, sample, 4)
            before = (Path(resumed) / "vq/batch-0000000.json").read_bytes()
            extend_batches(resumed, identity, 10, sample, 4)
            self.assertEqual(before, (Path(resumed) / "vq/batch-0000000.json").read_bytes())
            for arch in ("vq", "dima"):
                for index in range(3):
                    relative = f"{arch}/batch-{index:07d}.json"
                    self.assertEqual((Path(cold)/relative).read_bytes(),
                                     (Path(resumed)/relative).read_bytes())
            def no_recompute(*args):
                self.fail("cached pool resampled")
            self.assertEqual(extend_batches(resumed, identity, 10, no_recompute, 4),
                             {"vq": 10, "dima": 10})
            with self.assertRaises(ValueError):
                extend_batches(resumed, {"seed": 43}, 10, sample, 4)

if __name__ == "__main__":
    unittest.main()
