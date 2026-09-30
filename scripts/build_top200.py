"""Build the 200-candidate decision sheet from frozen outputs; no predictor rerun."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Indel
from amp_submission.io import AA, digest, read_fasta
from amp_submission.scoring import rewards
KD = dict(zip(AA, [1.8,2.5,-3.5,-3.5,2.8,-.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-.8,-.7,4.2,-.9,-1.3]))
OBJECTIVES = ("broad", "gram_positive", "gram_negative", "mdr", "selectivity", "joint")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("experiment", type=Path)
    ap.add_argument("--output", type=Path, default=Path("reports"))
    args=ap.parse_args()
    source=args.experiment/"challenge_20260929"
    out=args.output
    out.mkdir(exist_ok=True, parents=True)
    references=read_fasta(source/"official/data/antibacterial.fasta")
    cal=json.loads((source/"calibration.json").read_text())["scores"]
    rows=[]
    for arch in ("vq","dima"):
        folder=source/"deliverables"/arch
        top=read_fasta(folder/"generate/top.fasta")
        library=set(read_fasta(folder/"generate/library.fasta"))
        scored=[json.loads(x) for x in (folder/"top-scores.jsonl").read_text().splitlines()]
        assert len(top)==100 and [r["sequence"] for r in scored]==top
        nearest=np.load(source/f"sequence-audit/{arch}/top100/nearest_train_ratio.npy")
        assert read_fasta(source/f"sequence-audit/{arch}/top100/sequences.fasta")==top
        raw={k:np.array([r[k] for r in scored]) for k in cal}
        recalculated=rewards(raw,cal)
        for key in OBJECTIVES:
            assert np.allclose(recalculated[key],[r["reward_"+key] for r in scored])
        previous=[]
        for i, (sequence, data) in enumerate(zip(top,scored)):
            max_ref=process.extractOne(sequence,references,scorer=Indel.normalized_similarity)[1]
            max_prev=max((Indel.normalized_similarity(sequence,x) for x in previous), default=0.)
            row={"candidate_id":f"{arch.upper()}-{i+1:03d}", "architecture":arch,
                 "original_rank":i+1,"sequence":sequence,"length_aa":len(sequence),
                 "canonical_20aa":set(sequence)<=set(AA),"length_8_50":8<=len(sequence)<=50,
                 "member_of_own_library":sequence in library,"max_official_reference_ratio":max_ref,
                 "official_novelty_le_08":max_ref<=.8,"nearest_train_ratio":float(nearest[i]),
                 "train_exact_match":bool(nearest[i]>=1.),
                 "charge_proxy_KR_DE":sum(sequence.count(x) for x in "KR")-sum(sequence.count(x) for x in "DE"),
                 "hydrophobicity_Kyte_Doolittle":sum(KD[x] for x in sequence)/len(sequence),
                 **{k:v for k,v in data.items() if k!="sequence"},
                 "selection_utility_joint_minus_diversity":data["reward_joint"]-.15*max_prev,
                 "max_similarity_previously_selected":max_prev,
                 "linear_free_termini":"DESIGN_SPEC_ONLY_NOT_SYNTHESIS_VERIFIED",
                 "official_overall_success_rate":"NOT_MEASURED",
                 "official_gram_positive_success_rate":"NOT_MEASURED",
                 "official_gram_negative_success_rate":"NOT_MEASURED",
                 "official_mdr_success_rate":"NOT_MEASURED",
                 "official_MIC50_uM":"NOT_MEASURED","official_MIC90_uM":"NOT_MEASURED",
                 "official_HC50_uM":"NOT_MEASURED","official_safety_window":"NOT_MEASURED",
                 "official_any_MIC_le16uM":"NOT_MEASURED","synthesizability":"NOT_ASSESSED",
                 "independent_functional_oracle":"UNAVAILABLE",
                 "functional_coverage":"3 species proxies; 0/20 matched strains; 0/8 MDR-specific oracles"}
            rows.append(row);previous.append(sequence)
    for row in rows:
        row["unique_across_top200"]=sum(r["sequence"]==row["sequence"] for r in rows)==1
        row["nearest_other_top200_ratio"]=max(Indel.normalized_similarity(row["sequence"],r["sequence"])
                for r in rows if r["candidate_id"]!=row["candidate_id"])
    with (out/"top200-candidates.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/"top200-candidates.json").write_text(json.dumps(rows,indent=2)+"\n")
    summary={}
    for arch in ("vq","dima"):
        group=[r for r in rows if r["architecture"]==arch]
        metrics=["length_aa","max_official_reference_ratio","nearest_train_ratio",
                 "charge_proxy_KR_DE","hydrophobicity_Kyte_Doolittle",
                 "ania_ecoli","ania_paeruginosa","ania_saureus","hemopi2_hc50_um"]+["reward_"+x for x in OBJECTIVES]
        summary[arch]={key:{"min":float(min(r[key] for r in group)),
                           "median":float(np.median([r[key] for r in group])),
                           "mean":float(np.mean([r[key] for r in group])),
                           "max":float(max(r[key] for r in group))} for key in metrics}
        summary[arch]["checks"]={key:sum(bool(r[key]) for r in group) for key in
              ("canonical_20aa","length_8_50","member_of_own_library","official_novelty_le_08","train_exact_match","unique_across_top200")}
    report={"count":len(rows),"unique_count":len({r["sequence"] for r in rows}),
            "reference_sha256":digest(source/"official/data/antibacterial.fasta"),
            "calibration_sha256":digest(source/"calibration.json"),
            "summary":summary,
            "warning":"Proxy values are not official MIC, HC50 or competition results.",
            "sources":{arch:{n:digest(source/f"deliverables/{arch}"/n) for n in
                      ("generate/top.fasta","top-scores.jsonl")} for arch in ("vq","dima")}}
    (out/"top200-summary.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report))
if __name__=="__main__":
    main()
