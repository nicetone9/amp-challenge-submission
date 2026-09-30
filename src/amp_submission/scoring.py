"""External oracles retain their own licenses; no third-party weights are bundled."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path
import numpy as np

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
SPECIES = ("ecoli", "paeruginosa", "saureus")
ANIA_SUFFIX = {"ecoli": "EC", "paeruginosa": "PA", "saureus": "SA"}

def validate_sequences(sequences):
    if not sequences:
        raise ValueError("Empty candidate batch")
    for sequence in sequences:
        if not 8 <= len(sequence) <= 40:
            raise ValueError("Pilot requires length 8-40 so HemoPI2 does not truncate")
        if not set(sequence) <= VALID_AA:
            raise ValueError("Candidate contains non-canonical amino acids")


def predict_external(sequences, oracle_root, device):
    ANIA = Path(oracle_root).resolve() / 'ANIA'
    HEMO = Path(oracle_root).resolve() / 'HemoPI2'
    validate_sequences(sequences)
    ania_predict_code = (
        "import sys, torch; from torch.serialization import add_safe_globals; "
        "from torch.torch_version import TorchVersion; add_safe_globals([TorchVersion]); "
        "from src.inference.ania_predictor import predict_ania_from_csv; "
        "predict_ania_from_csv(sys.argv[1], sys.argv[2], sys.argv[3], device=sys.argv[4])"
    )

    ids = [f"seq_{i:06d}" for i in range(len(sequences))]
    id_to_index = {name: i for i, name in enumerate(ids)}
    result = {f"ania_{name}": np.full(len(sequences), np.nan) for name in SPECIES}
    result["hemopi2_hc50_um"] = np.full(len(sequences), np.nan)
    with tempfile.TemporaryDirectory(prefix="amp2-proxy-score-") as tmp_name:
        tmp = Path(tmp_name)
        fasta = tmp / "candidates.fasta"
        fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in zip(ids, sequences)), encoding="utf-8")
        cgr = tmp / "ania_features.csv"
        subprocess.run(
            [sys.executable, str(ANIA / "src/inference/fasta_encoder.py"), "--fasta", str(fasta), "--out", str(cgr)],
            cwd=ANIA, check=True, stdout=subprocess.DEVNULL,
        )
        for species, suffix in ANIA_SUFFIX.items():
            pred_path = tmp / f"ania_{suffix}.csv"
            subprocess.run(
                [
                    sys.executable, "-c", ania_predict_code, str(cgr),
                    str(ANIA / "weights" / f"ANIA_{suffix}.pt"), str(pred_path), str(device),
                ],
                cwd=ANIA, check=True, stdout=subprocess.DEVNULL,
            )
            with pred_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    result[f"ania_{species}"][id_to_index[row["ID"]]] = float(row["Predicted Log MIC Value"])

        hemo_run = subprocess.run(
            [sys.executable, str(HEMO / "code/hemopi2_regression.py"), "-i", str(fasta), "-wd", str(tmp)],
            cwd=HEMO, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        if hemo_run.returncode != 0:
            raise RuntimeError(f"HemoPI2 failed: {hemo_run.stdout[-1200:]}")
        hemo_path = tmp / "final_output.csv"
        with hemo_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                result["hemopi2_hc50_um"][id_to_index[row["SeqID"]]] = float(row["HC50(μM)"])

    if any(not np.isfinite(values).all() for values in result.values()):
        raise ValueError("External predictors did not return one finite result per sequence")
    return result

def rewards(scores, calibration):
    percentiles = {}
    for key, values in scores.items():
        c = np.sort(np.asarray(calibration[key]))
        x = np.asarray(values)
        p = (np.searchsorted(c,x,side="left")+np.searchsorted(c,x,side="right")) / (2*len(c))
        percentiles[key] = p if key == "hemopi2_hc50_um" else 1-p
    ec, pa, sa = [percentiles["ania_"+name] for name in SPECIES]
    broad = (ec+pa+sa)/3
    out = dict(broad=broad, gram_positive=sa, gram_negative=(ec+pa)/2,
               mdr=np.minimum.reduce([ec,pa,sa]),
               selectivity=.5*(broad+percentiles["hemopi2_hc50_um"]))
    out["joint"] = np.mean(list(out.values()),axis=0)
    return out
