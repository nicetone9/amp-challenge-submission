"""Exact cheap full-library statistics; expensive diversity metrics stay sample-labelled."""
import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path
import numpy as np
from amp_submission.io import AA, digest, read_fasta
from build_top200 import KD

def properties(seqs):
    lengths=np.array(list(map(len,seqs)))
    charge=np.array([sum(s.count(x) for x in "KR")-sum(s.count(x) for x in "DE") for s in seqs])
    hydro=np.array([sum(KD[x] for x in s)/len(s) for s in seqs])
    counts=Counter("".join(seqs));total=sum(counts.values())
    def summary(x):
        return dict(min=float(x.min()),q25=float(np.quantile(x,.25)),
                    median=float(np.median(x)),q75=float(np.quantile(x,.75)),
                    max=float(x.max()),mean=float(x.mean()))
    return {"n":len(seqs),"unique_n":len(set(seqs)),
            "legal_n":sum(8<=len(s)<=50 and set(s)<=set(AA) for s in seqs),
            "length":summary(lengths),"charge_KR_DE":summary(charge),"hydrophobicity_KD":summary(hydro),
            "length_histogram":dict(sorted(Counter(map(len,seqs)).items())),
            "aa_fraction":{a:counts[a]/total for a in AA},
            "repeat3_mean":float(np.mean([1-len({s[i:i+3] for i in range(len(s)-2)})/(len(s)-2) for s in seqs])),
            "repeat4_mean":float(np.mean([1-len({s[i:i+4] for i in range(len(s)-3)})/(len(s)-3) for s in seqs]))}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("experiment",type=Path)
    p.add_argument("--output",type=Path,default=Path("reports/library-statistics.json"))
    a=p.parse_args();source=a.experiment/"challenge_20260929"
    checks=json.loads((source/"submission-checks.json").read_text())
    train=set()
    with gzip.open(a.experiment/"prepared/manifest.csv.gz","rt") as handle:
        for row in csv.DictReader(handle):
            if row["split"]=="train":train.add(row["sequence"])
    libraries={arch:read_fasta(source/f"deliverables/{arch}/generate/library.fasta") for arch in ("vq","dima")}
    result={"scope":"all 50000 members per architecture; no subsampling for these statistics","groups":{}}
    for arch,seqs in libraries.items():
        path=source/f"deliverables/{arch}/generate/library.fasta"
        assert digest(path)==checks[arch]["verification"]["library_sha256"]
        prop=properties(seqs)
        prop["training_exact_hits"]=len(set(seqs)&train)
        prop["training_exact_hit_rate"]=prop["training_exact_hits"]/len(seqs)
        prop["official_reference_max_le_08"]="PASS from hash-matched existing full-library verification"
        prop["library_sha256"]=digest(path)
        prop["generation"]=checks[arch]["generation"]
        result["groups"][arch]=prop
    result["cross_architecture_shared_sequences"]=len(set(libraries["vq"])&set(libraries["dima"]))
    result["combined_unique_sequences"]=len(set(libraries["vq"])|set(libraries["dima"]))
    result["expensive_metrics_scope"]="CD50/CD95 and nearest-train/ProtT5 distributions remain sampled; see sequence-audit.json"
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result))
if __name__=="__main__":main()
