"""Publish aggregate audit evidence only; never upload protected AMP sequences."""
import argparse
import json
import os
from pathlib import Path
import swanlab
from .calibration import fingerprint
from .io import digest
from .prepare_reference import save_json

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("work/six-metrics"))
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--resume", action="store_true")
    a = p.parse_args()
    report_path = a.work / ("prepare-summary.json" if a.smoke else "calibration/report.json")
    report = json.loads(report_path.read_text())
    identity = {"report_sha256": digest(report_path), "pixi_lock_sha256": digest("pixi.lock"),
                "tracking_code_sha256": digest(__file__), "smoke": a.smoke,
                "scope": "aggregate screening audit, not model training or experimental activity"}
    run_id = fingerprint(identity)[:32]
    state_path = a.work / ("swanlab-smoke.json" if a.smoke else "swanlab-audit.json")
    previous = json.loads(state_path.read_text()) if state_path.exists() else None
    if previous and not a.resume:
        raise ValueError("Existing tracking state: use --resume")
    if a.resume and (not previous or previous["id"] != run_id):
        raise ValueError("Missing or mismatched tracking resume state")
    step = previous["last_step"] + 1 if previous else 0
    run = swanlab.init(workspace="nicetone9", project="AMP_step2challenge",
                       mode="online", id=run_id, resume="must" if a.resume else "allow",
                       name="six-metrics-" + ("smoke" if a.smoke else "partial-audit"),
                       group="six-metrics-seed42", job_type="reference-screening",
                       config=identity, log_dir=str(a.work / "swanlog"))
    metrics = {"audit/all_six_categories_passed": 0}
    if a.smoke:
        metrics.update({"preflight/ok": 1, "data/reference_eligible": report["reference_eligible"],
                        "data/candidate_count": report["pool_unique"],
                        "data/precheck_pass": report["pool_precheck_pass"]})
    else:
        for group in ("reference", "candidates"):
            for name, value in report[group]["per_metric_pass_rate"].items():
                metrics[group + "/pass_rate/" + name] = value
            metrics[group + "/joint_all_sequence_metrics"] = report[group]["joint_all_sequence_metrics"]
    swanlab.log(metrics, step=step)
    swanlab.log({"audit/status": swanlab.Text(report["status"] +
                 "; Missing metrics block completion. No wet-lab efficacy claim.")},
                step=1 if a.resume else 0)
    run.finish()
    save_json(a.work / ("swanlab-smoke.json" if a.smoke else "swanlab-audit.json"),
              {"id": run.id, "identity": identity, "last_step": step,
               "url": "https://swanlab.cn/@nicetone9/AMP_step2challenge/runs/" + run.id})
    print("SWANLAB_RUN_ID=" + run.id)

if __name__ == "__main__":
    main()
