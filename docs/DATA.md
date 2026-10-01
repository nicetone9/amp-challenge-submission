# Training data disclosure — exact prepared dataset PUBLIC

The complete prepared sequence dataset used by the two generator routes is now
[publicly downloadable](https://github.com/nicetone9/amp-challenge-submission/releases/download/frozen-pool-seed42-20260930/training-data-seed42-20260930.tar.gz).
[training-data.json](../training-data.json) records its size and SHA256.
Archive SHA256: `9ed5e6cf973831b806cebd02f5118287b4df6717103ae62488a702de2c22d70d`.
The original manifest is copied unchanged, not regenerated:
`a80a84884fa43c04ec95c082ddababf1297552e21d70d8f0a64d0693b2fddbbe`.
Frozen source milestone: AMP-Step0-v2-data-milestone-20260904, union input SHA256
`de30302dfd1d71d709c0d28f057cd578fb823b6bc629595af039d435a1a80a95`.

## Exact sequence and annotation scope

| Split | Core AMP | AMP-source candidate | Background / uncertain | Total |
| --- | ---: | ---: | ---: | ---: |
| Train | 9,212 | 525,292 | 531,601 | 1,066,105 |
| Validation | 1,202 | 65,607 | 66,112 | 132,921 |
| Internal generation evaluation | 1,185 | 64,563 | 67,990 | 133,738 |
| All | 11,599 | 655,462 | 665,703 | 1,332,764 |

All sequences are intact, canonical and 8–50 aa; no long-protein clipping.
Only `split=train` is base-model training data. Validation and internal
generation evaluation remain separately tagged; they are not protected downstream
ZYL test sets. All 50,252 rows with ZYL source membership are included.

The archive includes the exact manifest/audit, 4,963 antimicrobial and 2,159
antibacterial positive annotation rows actually joined onto that manifest
(task overlap is possible), original processing/training source snapshots,
source fingerprints, run protocols, licensing records and member checksums.
These are the annotation rows relevant to generator strata, not an export of
every raw task train CSV or unused target. Original scientific/provenance columns
are retained. Protected test/quarantine sequences are not published.

## Sources and processing

Common source screening preceded generator preparation: exact protection,
H30/cluster closure, noncanonical normalization and fixed motif exclusions.
The retained milestone record specifies MMseqs2 identity .30, coverage .60,
cov-mode 0, alignment-mode 3, seq-id-mode 2, sensitivity 8 and cluster-mode 2.
That historical screen was not rerun during release; its scope is the locked
reference/protocol, not a no-leakage guarantee for unknown test sets.

Preparation filters intact canonical 8–50 aa sequences, joins only eligible
positive train annotations, derives strata, shuffles sorted shared cluster IDs
with NumPy default_rng(42), assigns 80/10/10 clusters and stably sorts by
sequence_key. Source membership, component labels, OmegAMP curated headers,
cluster IDs and screen status are retained per row. sequence_key is an internal
sequence identity, not an original database accession.

Core: eligible ZYL antimicrobial/antibacterial positives or OmegAMP curated AMP.
Candidate: internally AMPSphere-named source membership excluding core.
Background is not an experimentally verified negative. Sampling is
core 30% / candidate 60% / background 10%, not row proportions or evidence grades.
ESM residue cache, train-only normalization and weighted train length distribution
are shared between architectures; ESM is original frozen ESM-2 650M, not Step-0
adapted weights.

The internally AMPSphere-named input is a four-CSV project package containing
UniProt background as well as AMP-related components, not a verified pure official
AMPSphere catalog. Original filenames and content hashes are disclosed; its
upstream accession/version map was not retained. Do not infer an official Zenodo
version from the internal name.

OmegAMP input hashes match the three public files at pinned commit
[c9500d3](https://github.com/szczurek-lab/OmegAMP/tree/c9500d3ca749f7db9f56409fc9868c2d09a62b3e/data/generative-model-data).
Its general component is marked Peptipedia/UniProt; the exact historical database
version and per-row Peptipedia-versus-UniProt map were not retained.

## Release verification

The archive's 29 members passed SHA256 verification. Public anonymous download
matches the archive hash. The original manifest checksum and split/stratum counts
match the training audit; reconstructing strata from the released positive rows
and curated-AMP component reproduces every stratum assignment. Prepared sequences
and shared clusters do not cross splits.

A new check against all 59 current protected ZYL test/quarantine files
(89,198 unique sequences) found zero exact overlap. Only their filenames,
fingerprints and aggregate counts are released, not their sequences.
This new exact check does not replace the historical H30/motif screening record.

## Rights, procedural data and remaining limits

The owner confirmed redistribution rights and authorized this release on
2026-09-30 (America/New_York). Owner-controlled annotation/curation rights are
offered under CC BY 4.0. Mixed prepared database rights retain ODbL 1.0 to preserve
potential Peptipedia share-alike obligations; individual third-party contents keep
their own terms. See [source-specific license scope](DATA_LICENSING.md).
This is not a blanket permissive relicense or a legal-clearance certificate.

RL policy updates use procedurally generated samples rather than another private
fixed training corpus. Original code, seeds, parent checkpoint hashes and
engineering-smoke/objective-proxy run reports are included. Exact historical
rollout strings/transitions were not saved and are not claimed recovered.
Externally pretrained ESM/ANIA/HemoPI2 upstream corpora are not repackaged as
project-owned data. ANIA/HemoPI2 overlap and calibration/generalization are not
independently established; ProtT5 is not an independent activity oracle.

The exact prepared dataset and used nonpublic annotations are now public.
Historical upstream accession/version gaps, procedural RL trajectory retention
and organizer interpretation of mixed data terms/fixed-pool reproduction remain
explicit limitations. Full co-authorship eligibility is not asserted.

Length and mass remain descriptive only for default selection. Legacy external
predictor routes retain their original 8–40 aa domain; a descriptor validity range
of 8–50 aa does not extend an external predictor's validated domain.
