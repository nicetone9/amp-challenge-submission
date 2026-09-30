"""Cold-install contracts; no model inference or external tools are invoked."""
import os
import runpy
import subprocess
import tomllib
import unittest
from pathlib import Path
from unittest.mock import patch

from amp_submission.generate import parser
from amp_submission.model_generate import run


class ColdInstallTests(unittest.TestCase):
    def test_r_is_only_a_legacy_dependency(self):
        project = tomllib.loads(Path("pyproject.toml").read_text())["project"]
        self.assertFalse(any(d.startswith("rpy2") for d in project["dependencies"]))
        self.assertIn("rpy2==3.5.17", project["optional-dependencies"]["legacy-oracles"])

    def test_missing_cdhit_fails_before_assets_or_sampling(self):
        with patch("amp_submission.model_generate.shutil.which", return_value=None), \
             patch("amp_submission.model_generate.verify_assets") as assets:
            with self.assertRaisesRegex(RuntimeError, "CD-HIT 4.8.1"):
                run(parser().parse_args(["--resample"]))
            assets.assert_not_called()

    def test_actual_cdhit_banner_passes_the_version_preflight(self):
        banner = "\t\t====== CD-HIT version 4.8.1 (built on Apr 24 2025) ======\n"
        result = subprocess.CompletedProcess(["cd-hit", "-h"], 1, banner)
        with patch("amp_submission.model_generate.shutil.which", return_value="/bin/cd-hit"), \
             patch("subprocess.run", return_value=result), \
             patch("amp_submission.model_generate.verify_assets", side_effect=RuntimeError("ASSETS_REACHED")):
            with self.assertRaisesRegex(RuntimeError, "ASSETS_REACHED"):
                run(parser().parse_args(["--resample"]))

    def test_cold_runtime_has_no_r_dependencies(self):
        config = tomllib.loads(Path("pixi.toml").read_text())
        self.assertEqual(set(config["feature"]["runtime"]["dependencies"]), {"python", "uv", "cd-hit"})
        self.assertTrue(config["environments"]["runtime"]["no-default-feature"])

    def test_cold_environment_does_not_inherit_preparing_venv(self):
        runner = runpy.run_path("scripts/reproduce_from_model.py")
        values = {"PATH": "/x/.pixi/envs/default/bin:/x/.venv/bin:/usr/bin",
                  "VIRTUAL_ENV": "/x/.venv", "PYTHONPATH": "/x/src", "R_HOME": "/x/R"}
        with patch.dict(os.environ, values, clear=True):
            env = runner["cold_environment"](Path("/fresh/run1"))
        self.assertEqual(env["PATH"], "/usr/bin")
        self.assertFalse({"VIRTUAL_ENV", "PYTHONPATH", "R_HOME"} & env.keys())
        self.assertEqual(env["UV_PROJECT_ENVIRONMENT"], "/fresh/run1/.venv")
        self.assertEqual(env["UV_CACHE_DIR"], "/fresh/run1/.uv-cache")

    def test_prepare_rejects_a_mutable_git_ref(self):
        runner = runpy.run_path("scripts/reproduce_from_model.py")
        with self.assertRaisesRegex(ValueError, "immutable Git commit SHA"):
            runner["prepare"](Path("/unused"), "main")

    def test_wrong_cdhit_version_fails_before_assets_or_sampling(self):
        result = subprocess.CompletedProcess(["cd-hit", "-h"], 0, "CD-HIT version 4.8.0\n")
        with patch("amp_submission.model_generate.shutil.which", return_value="/bin/cd-hit"), \
             patch("subprocess.run", return_value=result), \
             patch("amp_submission.model_generate.verify_assets") as assets:
            with self.assertRaisesRegex(RuntimeError, "CD-HIT 4.8.1"):
                run(parser().parse_args(["--resample"]))
            assets.assert_not_called()


if __name__ == "__main__":
    unittest.main()
