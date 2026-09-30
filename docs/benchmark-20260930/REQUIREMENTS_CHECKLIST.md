# Submission checklist — based on the requirements supplied by the team

This is a private preparation package, not a submitted or organizer-accepted entry.

## Minimum benchmark participation

| Requirement | Current status | Remaining action |
| --- | --- | --- |
| Method abstract | Prepared in ABSTRACT.md | Team review |
| 50,000 designed peptides | library.fasta: 50,000 unique, sequence-hard checks passed | Submit through the organizer's stated channel |
| Ranked Top100 and selection documentation | top.fasta: 100 from library; SELECTION_AND_MOTIFS.md and per-sequence CSVs | Team review of partial-quality ranking |
| Short data/resources/intervention summary | DATA_AND_INTERVENTIONS.md prepared | Confirm source records and disclose unresolved restrictions honestly |
| Private GitHub with model weights and inference code | Code is in nicetone9/amp-challenge-submission; weights exist on Lane but are not distributed | Resolve weight distribution and provide weights; verify clean inference |
| Read access for RasmusML and szymczakpau | Not granted or verified by this work | Repository owner grants read access and checks it |

The selected files and summary documents do not by themselves satisfy the missing weights/access requirement. Six-category quality completion is an internal evaluation goal; it is not represented here as an extra minimum requirement in the supplied text. Predictions are not wet-lab measurements.

## Full requirements / co-authorship eligibility

| Requirement | Current status | Remaining action |
| --- | --- | --- |
| Public template-based repository, weights, inference and detailed usage | Private repository and implementation exist; assets not distributed | Resolve release rights, publish only with authorization, complete cold-clone verification |
| Permissive OSI-approved license | MIT for original code/documentation | Preserve all separate upstream data/weight/code terms |
| Default seed and identical repeated generation | Seed 42; two fixed-pool selections match | Validate two full model-to-output runs via the competition entry point |
| Full training data disclosure and permissive release of nonpublic data | Partial disclosure, counts and fingerprints | Complete source/accession/version records and resolve/release permitted nonpublic data |
| Compliance with source terms | Not fully established | Finish source-by-source review before release |

The current legacy uv run generate entry point does not yet reproduce this mixed calibrated selection pipeline. Do not use the fixed-pool check as a substitute for full-entry-point reproducibility.

No collaborator permissions, repository visibility, data licenses, weight releases or competition submissions were changed while preparing these documents.
