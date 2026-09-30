# Training data, external resources and interventions

## Training input

This candidate package uses the frozen AMP-Step0-v2-data-milestone-20260904 data, not all subsequently updated Step 1 material. Existing audited disclosure records 1,332,764 eligible intact 8–50 aa canonical sequences: 1,066,105 train, 132,921 validation, and 133,738 internal generation evaluation. Core AMP sources are eligible antimicrobial/antibacterial positives from the collective dataset used in a previous study (internal ID: ZYL) and OmegAMP curated AMP; the candidate tier is AMPSphere excluding core; remaining eligible sequences are background/uncertain, not experimentally negative.

The shared cluster split uses seed 42 and 80/10/10 assignment, excluding protected test/quarantine sequences from the collective dataset (internal ID: ZYL). Training sampling is 30% core, 60% candidate, 10% background, not the raw data proportions. Both routes share original frozen ESM-2 650M residue representations and train-only normalization/length statistics.

For the current short-peptide calibration, 8,985 eligible training AMP sequences are used; sequence/embedding anchor and calibration-query clusters are separated. Motif discovery has its own fixed cluster discovery/validation split. Complete data accession/version/permission records and nonpublic-data release are still unresolved. See ../DATA.md for the existing counts and fingerprints. No private reference sequence list is included in this package.

## External resources

- Original facebook/esm2_t33_650M_UR50D: frozen representation model.
- DiMA upstream implementation: sequence-only latent diffusion components.
- ANIA: three species-level native log-MIC predictions used in policy rewards; no unsupported conversion to measured MIC in micromolar.
- HemoPI2: predicted HC50; not measured hemolysis.
- Official frozen antibacterial reference: 39,448 sequences for hard novelty audit.
- CD-HIT 4.8.1: deterministic clustering at identity 0.50, bidirectional coverage 0.80, single thread, retaining 8 aa peptides.
- ProtT5: planned independent embedding evaluation, not a completed criterion for this package.

Upstream training overlap, source terms, weight redistribution and oracle licensing require review. The original-code MIT license does not relicense external datasets or weights.

## Manual and computational interventions

No individual peptide was manually edited or substituted. Human-selected settings include the two architectures, predictor-reward objectives, frozen input milestone, seed, length range, low-complexity cutoff, clustering settings and high/low quotas. The user explicitly authorized lowering the internal quality percentile to obtain 50,000 candidates. The fixed grid was 0.50, 0.45, 0.40, 0.35, 0.30, 0.25; 0.30 was the highest feasible tested value. The initial pool contains equal raw counts from both generators; selection does not force the final source ratio to remain equal.

Reject noncanonical/out-of-range peptides, duplicates, single-residue fraction >0.6, exact known/training matches in the existing precheck, and whole-library official-reference similarity >0.8. Linear, unmodified, free termini are the design specification. No structural input or structural scoring is used.

This package's final ranking uses only the seven completed descriptor/novelty metrics, not the incomplete full MIC/safety score panel. Therefore policy-reward predictors and final candidate-ranking metrics must not be conflated.

## Final short-motif selection update

Ten Top100 slots (five high-tier and five low-tier) are reserved for frozen 3–4 aa motif matches; the remaining ninety preserve the same tier balance. The unchanged 50,000-member library passes the additional training-AMP novelty check. Safe Step 1 evidence is limited to 3,249 sequences intersecting the frozen AMP reference, with original holdout/quarantine flags excluded. This does not import all later Step 1 updates or establish coverage of every AMP database.
