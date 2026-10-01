# Training data disclosure — status PARTIAL

Frozen milestone: AMP-Step0-v2-data-milestone-20260904.
Sequence only, 20 canonical amino acids, intact length8–50, no long-protein clipping.
Milestone input SHA256: de30302dfd1d71d709c0d28f057cd578fb823b6bc629595af039d435a1a80a95.
Step2 manifest SHA256: a80a84884fa43c04ec95c082ddababf1297552e21d70d8f0a64d0693b2fddbbe.

| Split | Core AMP | AMP-source candidate | Background / uncertain | Total |
| --- | ---: | ---: | ---: | ---: |
| Train | 9,212 | 525,292 | 531,601 | 1,066,105 |
| Validation | 1,202 | 65,607 | 66,112 | 132,921 |
| Internal generation evaluation | 1,185 | 64,563 | 67,990 | 133,738 |
| All | 11,599 | 655,462 | 665,703 | 1,332,764 |

Sorted common clusters, NumPy default_rng42 shuffle,80/10/10 cluster assignment.
Audited sequence overlap0 and cluster overlap0. Protected test/quarantine sequences from the collective dataset (internal ID: ZYL)
excluded. Annotation joins only onto frozen milestone sequences; no Step1 resplit.
Sampling is core30%/candidate60%/background10%, not row proportions or evidence grades.

Core: eligible antimicrobial/antibacterial positives from the collective dataset
used in a previous study (internal ID: ZYL), or OmegAMP curated AMP.
Candidates: AMPSphere-source sequences excluding core.
Background is not interpreted as experimentally validated negative.
ESM residue cache, train-only normalization and weighted length distribution shared
between architectures. Default mixed-pool sampling/filtering uses the 8–50 aa validity range; length and
mass are descriptive, not quality gates. Explicit legacy ANIA/HemoPI2 predictor
routes retain their original 8–40 aa domain. Extending descriptor validity does
not extend an external predictor's validated domain.

## Missing for full disclosure

A name/count/hash is not a substitute for released training data. Exact accession
lists, source versions/URLs, complete processing records, each source's redistribution
terms and a permissively released nonpublic portion of the collective dataset remain to be cleared.
The repository does not currently satisfy the competition's full data-release condition.
No protected sequences are committed.

External ANIA/HemoPI2 training-data overlap and calibration/generalization are not
independently established. ProtT5 is an independent embedding model, not an independent
functional activity oracle.
