"""Frozen ProtT5 residue-mean embeddings, excluding EOS/PAD; shared reference/candidate path."""
import argparse
import json
import os
import time
from pathlib import Path
import numpy as np
import pandas as pd
from .calibration import fingerprint
from .io import digest
from .metric_scores import cosine_nearest
from .prepare_reference import save_json

def residue_mean(hidden, attention_mask):
    # T5 encodes one token per spaced standard amino acid followed by one EOS.
    import torch
    mask = attention_mask.bool().clone()
    eos = mask.sum(1) - 1
    mask[torch.arange(len(mask), device=mask.device), eos] = False
    if torch.any(mask.sum(1) == 0):
        raise ValueError("Empty residue mask")
    return (hidden.float() * mask[..., None]).sum(1) / mask.sum(1)[:, None]

def run(work, model_root, device="cuda", smoke=False, track=False):
    import torch
    from transformers import T5EncoderModel, T5Tokenizer
    work, model_root = Path(work), Path(model_root)
    info = json.loads((model_root / "fingerprint.json").read_text())
    for name, expected in info["sha256"].items():
        if digest(model_root / name) != expected:
            raise ValueError("Changed ProtT5 asset: " + name)
    prepared = json.loads((work / "prepare-summary.json").read_text())
    for name in ("reference.csv", "candidates.csv"):
        if digest(work / name) != prepared["artifacts"][name]:
            raise ValueError("Changed embedding input table")
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA allocation required; no CPU fallback")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    torch.set_num_threads(4)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.benchmark = False
    identity = {"model": info, "code_sha256": digest(__file__), "device": device,
                "reference_identity": digest(work / "identity.json"),
                "pixi_lock": digest("pixi.lock"), "batch_size": 32, "smoke": smoke,
                "dtype": "float16" if device == "cuda" else "float32"}
    target = work / ("embedding-smoke" if smoke else "scores/embedding")
    if (target / "identity.json").exists() and json.loads((target / "identity.json").read_text()) != identity:
        raise ValueError("Embedding resume fingerprint mismatch")
    save_json(target / "identity.json", identity)
    run_handle, tracking = None, None
    if track:
        import swanlab
        tracking_path = target / "tracking.json"
        tracking = json.loads(tracking_path.read_text()) if tracking_path.exists() else None
        run_id = fingerprint([identity, str(target.resolve())])[:32]
        if tracking and tracking["id"] != run_id:
            raise ValueError("Embedding tracking resume identity mismatch")
        run_handle = swanlab.init(workspace="nicetone9", project="AMP_step2challenge",
            mode="online", id=run_id, resume="must" if tracking else "allow",
            name="six-metrics-prott5", group="six-metrics-seed42", job_type="embedding",
            config={"identity_sha256": fingerprint(identity), "device": device,
                    "batch_size": 32, "scope": "independent sequence embedding evaluation"},
            log_dir=str(work / "swanlog"))
        tracking = tracking or {"id": run_id, "last_step": -1}
        save_json(tracking_path, tracking)
    tokenizer = T5Tokenizer.from_pretrained(model_root, local_files_only=True, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_root, local_files_only=True,
                torch_dtype=torch.float16 if device == "cuda" else torch.float32).to(device).eval()
    model.requires_grad_(False)
    started = time.monotonic()
    tables, vectors = {}, {}
    for name in ("reference", "candidates"):
        table = pd.read_csv(work / (name + ".csv"), keep_default_na=False)
        if smoke:
            table = table.iloc[:2].copy()
        tables[name] = table
        rows = []
        for start in range(0, len(table), 32):
            sequences = table.sequence.iloc[start:start+32].tolist()
            key = fingerprint({"identity": identity, "sequences": sequences})
            path = target / "cache" / (key + ".npz")
            meta = path.with_suffix(".json")
            if meta.exists():
                cached = json.loads(meta.read_text())
                if cached["identity"] != key or digest(path) != cached["sha256"]:
                    raise ValueError("Corrupt embedding batch cache")
                values = np.load(path, allow_pickle=False)["embeddings"]
            else:
                tokens = tokenizer([" ".join(s) for s in sequences], return_tensors="pt", padding=True)
                if tokens.attention_mask.sum(1).tolist() != [len(s)+1 for s in sequences]:
                    raise ValueError("Tokenizer/residue alignment mismatch")
                tokens = {k: v.to(device) for k, v in tokens.items()}
                with torch.inference_mode():
                    hidden = model(**tokens).last_hidden_state
                    values = residue_mean(hidden, tokens["attention_mask"]).cpu().numpy()
                path.parent.mkdir(parents=True, exist_ok=True)
                np.savez(path, embeddings=values)
                save_json(meta, {"identity": key, "sha256": digest(path)})
            if values.shape != (len(sequences), 1024) or not np.isfinite(values).all():
                raise ValueError("Invalid ProtT5 embedding output")
            rows.append(values)
            count = start + len(sequences)
            if track and (count % 1024 == 0 or count == len(table)):
                tracking["last_step"] += 1
                swanlab.log({"embedding/" + name + "_completed": count,
                             "embedding/wall_seconds": time.monotonic() - started},
                            step=tracking["last_step"])
                save_json(tracking_path, tracking)
        vectors[name] = np.concatenate(rows)
        np.savez(target / (name + "-vectors.npz"), embeddings=vectors[name],
                 sequences=np.asarray(table.sequence.tolist()))
        print(json.dumps({"table": name, "embedded": len(table), "smoke": smoke}), flush=True)
    if smoke:
        report = {"status": "SMOKE_ONLY", "reference_n": len(tables["reference"]),
                  "candidate_n": len(tables["candidates"]), "identity": fingerprint(identity)}
    else:
        roles = tables["reference"].reference_role
        anchors = vectors["reference"][roles.eq("anchor")]
        query = roles.eq("calibration")
        for name in tables:
            indices = np.flatnonzero(query) if name == "reference" else np.arange(len(tables[name]))
            values = cosine_nearest(vectors[name][indices], anchors)
            pd.DataFrame({"sequence": tables[name].sequence.iloc[indices].tolist(),
                          "embedding_nn_cosine_distance": values}).to_csv(target / (name + ".csv"), index=False)
        report = {"status": "COMPUTED", "identity": fingerprint(identity),
                  "sha256": {name: digest(target / (name + ".csv")) for name in tables},
                  "vectors_sha256": {name: digest(target / (name + "-vectors.npz")) for name in tables}}
    report["wall_seconds"] = time.monotonic() - started
    save_json(target / "complete.json", report)
    if run_handle is not None:
        run_handle.finish()
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, default=Path("work/six-metrics"))
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--track", action="store_true")
    a = p.parse_args()
    print(json.dumps(run(a.work, a.model, a.device, a.smoke, a.track)))

if __name__ == "__main__":
    main()
