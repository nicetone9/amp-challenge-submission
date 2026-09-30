"""Resumable shared scoring; batch caches bind ordered sequences, tools and reference."""
import argparse
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
from .calibration import REGISTRY, calibrate, fingerprint, registry_payload
from .io import digest
from .metric_scores import descriptors, anchor_novelty
from .prepare_reference import save_json
from .scoring import predict_external

def cache_scores(sequences, directory, identity, compute, batch_size=512, progress_callback=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for start in range(0, len(sequences), batch_size):
        batch = sequences[start:start+batch_size]
        key = fingerprint({"identity": identity, "sequences": batch, "batch_size": batch_size})
        path = directory / (key + ".json")
        cache_hit = path.exists()
        if cache_hit:
            cached = json.loads(path.read_text())
            if cached["key"] != key or cached["score_sha256"] != fingerprint(cached["scores"]):
                raise ValueError("Corrupt score cache")
            scores = cached["scores"]
        else:
            values = compute(batch)
            scores = {name: np.asarray(value, float).tolist() for name, value in values.items()}
            if any(len(value) != len(batch) or not np.isfinite(value).all() for value in scores.values()):
                raise ValueError("Scorer returned missing, nonfinite or misaligned values")
            save_json(path, {"key": key, "scores": scores, "score_sha256": fingerprint(scores)})
        if outputs and set(scores) != set(outputs):
            raise ValueError("Inconsistent score columns")
        for name, values in scores.items():
            outputs.setdefault(name, []).extend(values)
        if progress_callback is not None:
            progress_callback(start + len(batch), cache_hit)
    return {name: np.asarray(values) for name, values in outputs.items()}

def oracle_fingerprint(root):
    root = Path(root)
    files = []
    for relative in ("ANIA/src", "ANIA/configs", "ANIA/weights", "HemoPI2/code", "HemoPI2/Model"):
        folder = root / relative
        if not folder.is_dir():
            raise ValueError("Missing oracle directory: " + relative)
        files.extend(p for p in folder.rglob("*") if p.is_file()
                     and "__pycache__" not in p.parts and p.suffix != ".pyc")
    return {str(p.relative_to(root)): digest(p) for p in sorted(files)}

def score_stage(work, stage, oracles=None, device="cpu", track=False):
    work = Path(work)
    prepared = json.loads((work / "prepare-summary.json").read_text())
    for name, expected in prepared["artifacts"].items():
        if digest(work / name) != expected:
            raise ValueError("Prepared artifact changed: " + name)
    tables = {name: pd.read_csv(work / (name + ".csv"), keep_default_na=False)
              for name in ("reference", "candidates")}
    anchor = tables["reference"].query("reference_role == 'anchor'").sequence.tolist()
    tools = oracle_fingerprint(oracles) if stage == "oracles" else {}
    identity = {"reference_identity": digest(work / "identity.json"),
                "registry": fingerprint(registry_payload()), "stage": stage, "device": device,
                "code": {name: digest(Path(__file__).parent / name) for name in
                         ("score_stage.py", "metric_scores.py", "scoring.py", "calibration.py")},
                "pixi_lock": digest("pixi.lock"), "tool_assets": tools, "batch_size": 512}
    target = work / "scores" / stage
    previous = target / "identity.json"
    if previous.exists() and json.loads(previous.read_text()) != identity:
        raise ValueError("Changed scoring environment: use a new work directory")
    save_json(previous, identity)
    started = time.monotonic()
    run, tracking = None, None
    if track:
        import swanlab
        tracking_path = target / "tracking.json"
        tracking = json.loads(tracking_path.read_text()) if tracking_path.exists() else None
        run_id = fingerprint([identity, str(work.resolve())])[:32]
        if tracking and tracking["id"] != run_id:
            raise ValueError("Tracking resume identity mismatch")
        run = swanlab.init(workspace="nicetone9", project="AMP_step2challenge", mode="online",
                          id=run_id, resume="must" if tracking else "allow",
                          name="six-metrics-" + stage, group="six-metrics-seed42",
                          job_type="full-scoring", config={"identity_sha256": fingerprint(identity),
                          "stage": stage, "device": device, "seed": 42, "batch_size": 512},
                          log_dir=str(work / "swanlog"))
        tracking = tracking or {"id": run_id, "last_step": -1}
        save_json(tracking_path, tracking)
    for name, table in tables.items():
        if stage == "novelty" and name == "reference":
            table = table[table.reference_role.eq("calibration")].copy()
        sequences = table.sequence.tolist()
        if stage == "descriptors":
            compute = descriptors
        elif stage == "novelty":
            compute = lambda seqs: {"anchor_novelty": anchor_novelty(seqs, anchor)}
        else:
            compute = lambda seqs: predict_external(seqs, oracles, device)
        def progress(count, cache_hit):
            if run is None or cache_hit or (count % 2560 and count != len(sequences)):
                return
            tracking["last_step"] += 1
            swanlab.log({"scoring/" + name + "_completed": count,
                         "scoring/wall_seconds": time.monotonic() - started},
                        step=tracking["last_step"])
            save_json(tracking_path, tracking)
        try:
            scores = cache_scores(sequences, work / "cache" / stage, identity, compute,
                                  progress_callback=progress)
        except BaseException as error:
            if run is not None:
                run.finish(state="crashed", error=type(error).__name__ + ": " + str(error))
            raise
        result = pd.DataFrame({"sequence": sequences, **scores})
        result.to_csv(target / (name + ".csv"), index=False)
        print(json.dumps({"stage": stage, "table": name, "rows": len(result)}), flush=True)
    report = {"status": "COMPUTED", "stage": stage, "identity": fingerprint(identity),
              "wall_seconds": time.monotonic() - started,
              "sha256": {name: digest(target / (name + ".csv")) for name in tables}}
    save_json(target / "complete.json", report)
    if run is not None:
        run.finish()
    return report

def load_scores(work, name):
    frame = pd.read_csv(work / (name + ".csv"), keep_default_na=False)
    source_hashes = {}
    for folder in sorted((work / "scores").glob("*")):
        if not (folder / "complete.json").exists():
            continue
        completion = json.loads((folder / "complete.json").read_text())
        path = folder / (name + ".csv")
        if digest(path) != completion["sha256"][name]:
            raise ValueError("Scoring artifact hash mismatch")
        scores = pd.read_csv(path)
        if not scores.sequence.is_unique:
            raise ValueError("Non-unique score rows")
        overlap = set(scores.columns[1:]) & set(frame.columns)
        if overlap == {"length"}:
            expected = scores.sequence.str.len().to_numpy()
            if not np.array_equal(scores["length"].to_numpy(), expected):
                raise ValueError("Scored lengths differ from full sequence lengths")
            frame = frame.drop(columns=["length"])
        elif overlap:
            raise ValueError("Duplicate metric column")
        frame = frame.merge(scores, on="sequence", how="left", validate="one_to_one", sort=False)
        source_hashes[folder.name] = digest(folder / "complete.json")
    return frame, source_hashes

def audit_calibration(work, threshold=.5):
    work = Path(work)
    ref, r_hash = load_scores(work, "reference")
    cand, c_hash = load_scores(work, "candidates")
    ref = ref[ref.reference_role.eq("calibration")].copy()
    metric_names = [m.name for m in REGISTRY if m.level == "sequence"]
    r_scores = {name: ref[name].to_numpy(float) for name in metric_names if name in ref}
    c_scores = {name: cand[name].to_numpy(float) for name in metric_names if name in cand}
    reference = calibrate(r_scores, r_scores, len(ref), threshold)
    candidate = calibrate(r_scores, c_scores, len(cand), threshold)
    frozen = {"reference_sha256": digest(work / "reference.csv"),
              "registry": registry_payload(), "threshold": threshold,
              "reference_score_sources": r_hash,
              "scores": {name: values.tolist() for name, values in r_scores.items()}}
    caldir = work / "calibration"
    ready = all(value == "COMPUTED" for value in candidate["statuses"].values())
    path = caldir / ("frozen.json" if ready else "partial-reference.json")
    if ready and path.exists() and json.loads(path.read_text()) != frozen:
        raise ValueError("Frozen reference changed: choose a new calibration directory")
    save_json(path, frozen)
    for name, p in candidate["percentiles"].items():
        cand[name + "_percentile"] = p
        cand[name + "_pass"] = np.isfinite(p) & (p >= threshold)
    cand["all_sequence_metrics_pass"] = candidate["passed"]
    cand["internal_ranking_score"] = candidate["ranking"]
    cand["weakest_percentile"] = candidate["weakest"]
    # This is NOT submission readiness: official-reference and set gates are separate.
    cand["qualified_for_next_audit"] = candidate["passed"] & cand.hard_precheck.astype(bool)
    cand.to_csv(caldir / "candidate-percentiles.csv", index=False)
    available = [name for name, state in candidate["statuses"].items() if state == "COMPUTED"]
    def summary(out):
        p = out["percentiles"]
        return {"per_metric_pass_rate": {name: float((np.isfinite(values) & (values >= threshold)).mean())
                                         for name, values in p.items()},
                "joint_all_sequence_metrics": float(out["passed"].mean()),
                "joint_computed_metrics_only_DIAGNOSTIC": float(np.logical_and.reduce(
                    [np.isfinite(p[name]) & (p[name] >= threshold) for name in available]).mean())
                    if available else None}
    report = {"status": "INCOMPLETE", "threshold": threshold, "reference_n": len(ref),
              "candidate_n": len(cand), "reference": summary(reference), "candidates": summary(candidate),
              "metric_status": candidate["statuses"], "frozen_sha256": digest(path),
              "score_sources": c_hash, "all_six_categories_passed": False,
              "pending": ["official-reference hard audit", "collection-level gates",
                          "motif discovery/validation", "deterministic full generation"],
              "note": "Computed-only diagnostic cannot qualify candidates; missing metrics fail closed."}
    save_json(caldir / "report.json", report)
    print(json.dumps(report), flush=True)
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("descriptors", "novelty", "oracles", "calibrate"))
    p.add_argument("--work", type=Path, default=Path("work/six-metrics"))
    p.add_argument("--oracles", type=Path, default=Path("external"))
    p.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    p.add_argument("--reference-percentile", type=float, default=.5)
    p.add_argument("--track", action="store_true")
    a = p.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    if a.stage == "calibrate":
        audit_calibration(a.work, a.reference_percentile)
    else:
        score_stage(a.work, a.stage, a.oracles, a.device, a.track)

if __name__ == "__main__":
    main()
