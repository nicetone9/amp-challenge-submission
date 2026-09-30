"""Track only verified aggregate evidence for the length-hard-only delivery."""
import json
from pathlib import Path
import swanlab
from amp_submission.io import digest
from amp_submission.prepare_reference import save_json

report_path = Path("reports/length-hard-only-seed42-verification.json")
report = json.loads(report_path.read_text())
if report["status"] != "LENGTH_HARD_ONLY_FIXED_POOL_DOUBLE_RUN_VERIFIED":
    raise ValueError("Selection verification has not passed")
identity = {"report_sha256": digest(report_path), "seed": report["seed"],
            "quality_metrics": report["quality_metrics"],
            "descriptive_only": report["descriptive_only"],
            "fasta_sha256": report["fasta_sha256"], "inputs": report["inputs"],
            "scope": report["scope"], "length_hard_bounds": [8, 50]}
started = Path("work/length-hard-only-seed42/swanlab-started.json")
if started.exists():
    raise ValueError("Tracking already attempted; inspect state and use an explicit resume workflow")
save_json(started, identity)
run = swanlab.init(workspace="nicetone9", project="AMP_step2challenge",
                   mode="online", id=identity["report_sha256"][:32], resume="allow",
                   name="length-hard-only-fixed-pool", group="length-hard-only-seed42",
                   job_type="selection-audit", config=identity,
                   log_dir="work/length-hard-only-seed42/swanlog")
metrics = {"audit/byte_identical": 1, "audit/known_reference_hard_pass": 1,
           "audit/all_six_categories_passed": 0,
           "selection/reference_percentile": report["selected_threshold"]}
for name, stats in report["statistics"].items():
    for key in ("count", "min", "max", "mean", "median", "clusters", "motif_supported"):
        metrics[name + "/" + key] = stats[key]
    for key, value in stats["bins"].items():
        metrics[name + "/length/" + key] = value
swanlab.log(metrics, step=0)
swanlab.log({"audit/status": swanlab.Text(report["scope"] + "; " + report["warning"])}, step=0)
run.finish()
save_json("work/length-hard-only-seed42/swanlab.json",
          {"run_id": run.id, "report_sha256": identity["report_sha256"],
           "url": "https://swanlab.cn/@nicetone9/AMP_step2challenge/runs/" + run.id})
print("SWANLAB_RUN_ID=" + run.id, flush=True)
