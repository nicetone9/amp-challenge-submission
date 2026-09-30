"""Run the unmodified default command twice in isolated cold candidate workspaces."""
import argparse
import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import time

from amp_submission.io import digest, read_fasta, verify_assets
from amp_submission.prepare_reference import save_json


REPOSITORY = "https://github.com/nicetone9/amp-challenge-submission.git"


def cold_environment(target):
    env = dict(os.environ, UV_CACHE_DIR=str(target / ".uv-cache"),
               UV_PROJECT_ENVIRONMENT=str(target / ".venv"), PYTHONNOUSERSITE="1")
    for key in ("VIRTUAL_ENV", "PYTHONPATH", "PYTHONHOME", "R_HOME"):
        env.pop(key, None)
    # Do not inherit the preparing checkout's Python/R/tool environment.
    env["PATH"] = os.pathsep.join(p for p in env["PATH"].split(os.pathsep)
                                 if "/.pixi/envs/" not in p and "/.venv/" not in p)
    return env


def runtime_command(target):
    return ["pixi", "run", "--locked", "-e", "runtime",
            "--manifest-path", str(target / "pixi.toml")]


def prepare(root, revision):
    source = Path.cwd().resolve()
    if not revision or len(revision) != 40 or set(revision) - set("0123456789abcdef"):
        raise ValueError("Preparation requires a full immutable Git commit SHA")
    if root.exists():
        raise ValueError("Use a new reproducibility directory")
    bundles = {relative: verify_assets(source / relative) for relative in
               ("checkpoint", "checkpoint/generation-reference")}
    paths = ["pyproject.toml", "uv.lock", "pixi.toml", "pixi.lock", ".python-version"]
    paths += [str(p) for p in sorted(Path("src").rglob("*.py"))]
    paths += ["scripts/verify_submission.py", "scripts/reproduce_from_model.py"]
    hashes = {name: digest(name) for name in paths}
    root.mkdir(parents=True)
    installations = {}
    for name in ("run1", "run2"):
        target = root / name
        env = cold_environment(target)
        with (root / (name + "-install.log")).open("w") as log:
            subprocess.run(["git", "clone", "--no-checkout", REPOSITORY, str(target)],
                           check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
            subprocess.run(["git", "checkout", "--detach", revision], cwd=target,
                           check=True, stdout=log, stderr=subprocess.STDOUT, env=env)
            actual = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=target,
                                             text=True, env=env).strip()
            if actual != revision or any(digest(target / p) != h for p, h in hashes.items()):
                raise ValueError("Cloned revision differs from the audited source/locks")
            # Provision only manifest-listed frozen inputs, never generated candidates.
            # Copies are isolated; there is no symlink to the staging checkout.
            for relative, manifest in bundles.items():
                for filename in ("manifest.json", *manifest["sha256"]):
                    dst = target / relative / filename
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(source / relative / filename, dst)
            subprocess.run([*runtime_command(target), "uv", "sync", "--frozen"],
                           cwd=target, check=True, env=env, stdout=log, stderr=subprocess.STDOUT)
            probe_code = (
                "import importlib.util,json,shutil,sys; import torch; "
                "print(json.dumps(dict(python=sys.version,prefix=sys.prefix,"
                "rpy2_installed=importlib.util.find_spec('rpy2') is not None,"
                "r_on_path=shutil.which('R'),cd_hit=shutil.which('cd-hit'),torch=torch.__version__)))"
            )
            probe = subprocess.run([*runtime_command(target), "uv", "run", "--frozen",
                                    "python", "-c", probe_code], cwd=target, check=True,
                                   env=env, stdout=subprocess.PIPE, stderr=log, text=True)
            installations[name] = json.loads(probe.stdout)
            if (installations[name]["rpy2_installed"] or
                    Path(installations[name]["prefix"]).resolve() != target / ".venv"):
                raise ValueError("Cold default installation inherited legacy Python dependencies")
        if (target / "generate").exists() or list(target.glob(".amp-fresh-*")):
            raise ValueError("Preparation unexpectedly contains generated output")
    report = {"status": "PREPARED_NOT_EXECUTED", "source": str(source), "source_files": hashes,
              "repository": REPOSITORY, "revision": revision, "installations": installations,
              "asset_delivery": "private manifest-listed copies; public distribution unresolved",
              "assets_manifest": digest(source / "checkpoint/manifest.json"),
              "reference_manifest": digest(source / "checkpoint/generation-reference/manifest.json"),
              "command": ["uv", "run", "generate"],
              "native_prerequisites": "fresh locked Pixi runtime per clone; Python, uv, CD-HIT"}
    save_json(root / "prepared.json", report)
    print(json.dumps(report), flush=True)


def check_inputs(root, prepared):
    for name in ("run1", "run2"):
        target = root / name
        if prepared.get("revision") and subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=target, text=True).strip() != prepared["revision"]:
            raise ValueError("Changed cloned revision")
        for relative, expected in prepared["source_files"].items():
            if digest(target / relative) != expected:
                raise ValueError("Changed isolated source: " + name + "/" + relative)
        for relative, key in (("checkpoint", "assets_manifest"),
                              ("checkpoint/generation-reference", "reference_manifest")):
            verify_assets(target / relative)
            if digest(target / relative / "manifest.json") != prepared[key]:
                raise ValueError("Changed frozen assets")


def verify(root, prepared):
    reports = [json.loads((root / name / "generate/generation.json").read_text())
               for name in ("run1", "run2")]
    for key in ("seed", "pool_size", "batch_size", "threads", "n_sequences", "top_k",
                "motif_quota", "thresholds", "assets_sha256", "reference_sha256",
                "code_sha256", "locks", "raw_sha256", "environment"):
        if reports[0][key] != reports[1][key]:
            raise ValueError("Two-run input/raw/environment mismatch: " + key)
    if reports[0]["artifact_directory"] == reports[1]["artifact_directory"]:
        raise ValueError("Two runs shared candidate work")
    for report in reports:
        if not report["fresh_model_sampling"] or report["candidate_cache_read"]:
            raise ValueError("Not a cold model generation")
        if (report["n_sequences"], report["top_k"]) != (50000, 100):
            raise ValueError("Nondefault quantities cannot certify the submission")
    artifacts = {}
    for name in ("library.fasta", "top.fasta"):
        first, second = (root / run / "generate" / name for run in ("run1", "run2"))
        if first.read_bytes() != second.read_bytes():
            raise ValueError("Byte-for-byte reproducibility failed: " + name)
        artifacts[name] = digest(first)
    validators = runpy.run_path(str(root / "run1/scripts/verify_submission.py"))
    references = set(read_fasta(root / "run1/checkpoint/antibacterial.fasta"))
    for run in ("run1", "run2"):
        output = root / run / "generate"
        library = validators["_verify_sequences"](output / "library.fasta")
        validators["_verify_top"](output / "top.fasta", library, 100)
        validators["_verify_no_overlap"](library, references)
        validators["_veritfy_max_simularity"](set(read_fasta(output / "top.fasta")), references)
    check_inputs(root, prepared)
    result = {"status": "TWO_COLD_MODEL_RUNS_BYTE_IDENTICAL",
              "command_each_run": ["uv", "run", "generate"], "seed": 42,
              "n_sequences": 50000, "top_k": 100, "fasta_sha256": artifacts,
              "raw_sha256": reports[0]["raw_sha256"], "environment": reports[0]["environment"],
              "inputs": prepared, "run_reports_sha256": [
                  digest(root / name / "generate/generation.json") for name in ("run1", "run2")],
              "official_template_component_checks_passed": True,
              "swanlab_run_ids": [r["swanlab_run_id"] for r in reports],
              "not_claimed": ["cross-hardware identity", "organizer acceptance",
                              "public asset delivery or complete official cold-clone acceptance",
                              "six-category quality completion", "wet-lab activity"]}
    save_json(root / "verification.json", result)
    print(json.dumps(result), flush=True)
    return result


def execute(root):
    prepared = json.loads((root / "prepared.json").read_text())
    check_inputs(root, prepared)
    settings = dict(UV_FROZEN="1", PYTHONHASHSEED="42",
               OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4",
               CUBLAS_WORKSPACE_CONFIG=":4096:8", TOKENIZERS_PARALLELISM="false",
               AMP_TRACK_ONLINE="1")
    # Fail before either invocation if the run directories are not cold.
    for name in ("run1", "run2"):
        if ((root / name / "generate").exists() or
                list((root / name).glob(".amp-fresh-*"))):
            raise ValueError("Cold run directory already used: " + name)
    for name in ("run1", "run2"):
        env = dict(cold_environment(root / name), **settings)
        start = time.monotonic()
        save_json(root / "status.json", {"status": "RUNNING", "run": name})
        with (root / (name + ".log")).open("w") as log:
            result = subprocess.run([*runtime_command(root / name), "uv", "run", "generate"], cwd=root / name,
                                    env=env, stdout=log, stderr=subprocess.STDOUT)
        save_json(root / (name + "-exit.json"),
                  {"exit_code": result.returncode, "wall_seconds": time.monotonic()-start})
        if result.returncode:
            save_json(root / "status.json", {"status": "FAILED", "run": name,
                                             "exit_code": result.returncode})
            raise RuntimeError(name + " failed; dependent run was not started")
    verify(root, prepared)
    save_json(root / "status.json", {"status": "VERIFIED"})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("work/from-model-repro-seed42"))
    p.add_argument("--phase", choices=("prepare", "run"), required=True)
    p.add_argument("--revision", help="Full Git commit SHA; required for preparation")
    a = p.parse_args()
    root = a.root.resolve()
    if a.phase == "prepare":
        prepare(root, a.revision)
    else:
        execute(root)
