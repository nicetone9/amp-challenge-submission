# Competition alignment (checked 2026-09-30)

Sources: [official template](https://github.com/szczurek-lab/amp-challenge-2027),
[category definitions](https://szczurek-lab.github.io/amp-challenge-website/#categories),
[Kaggle overview](https://www.kaggle.com/competitions/amp-challenge/overview).

## Five experimental categories

| Official category | Primary criterion | Tie-break / gate | Current evidence |
| --- | --- | --- | --- |
| Broad-spectrum | Success fraction across20strains | MIC90 over20 | EC,PA,SA species-level proxies only; no20-strain MIC panel |
| Gram-positive | Success fraction across5Gram+strains | MIC50 overGram+ | SA proxy only; no5-strain measurements |
| Gram-negative | Success fraction across15Gram−strains | MIC50 overGram− | EC,PA proxies only; no15-strain measurements |
| MDR ESKAPE | Success fraction across8MDRisolates | MIC50 overMDR | NoMDR-specific oracle; worst-of-three proxy only |
| Optimal selectivity | HC50/MIC50 across20strains | AtleastoneMIC≤16µM; inactiveexcluded; cappedHC50tie broken byMIC50 | HC50 prediction+rank proxy only; official safetywindow unavailable |

Activity means measured MIC≤16µM. MIC assay limit64µM; no inhibition is>64.
HC50 assay limit128µM; nonhemolytic is>128,not128.
Team results average peptide metrics over25randomly sampled candidates fromTop100.
Individual-peptide proxy scores cannot establish team category results.

## Current sequence and delivery gates

The current delivery is the new seed-42 300,000 VQ + 300,000 DiMA pool, not the
historical separate-architecture libraries. Minimum and full requirements follow
the [official template](https://github.com/szczurek-lab/amp-challenge-2027#submission-requirements).

| Requirement | Current evidence / remaining gate |
| --- | --- |
| Method abstract and ranking/filter summary | Current ABSTRACT.md and WRITEUP.md |
| Library of 50,000 unique designs | Source finalization and independent metric table pass |
| Ranked Top100 from that library | 100 unique members; 50 per tier and ten motif quota selections |
| 20 standard AA; length 8–50 | Library 8–50; Top100 8–36; checked |
| Linear, free termini, unmodified | Design specification; synthesis/QC not assessed |
| No exact official-reference library overlap | Full-library combined-reference validation passes |
| Top100 reference Levenshtein ratio <=0.8 | Stricter gate also applied to the entire library |
| GitHub inference code and weights | Public code and approved frozen inference release |
| Organizer read access | Public repository/assets can be read anonymously; no collaborator/admin action required |
| Permissive original-code license | MIT; DiMA/ESM MIT and validator BSD notices retained separately |
| Defined Python and uv.lock | Python 3.11; locked install and default entry implemented |
| No-argument uv run generate | Fixed-pool filtering; native CD-HIT prerequisite is supplied by locked Pixi runtime |
| Identical repeated default output | [Two public cold clones pass](../reports/frozen-pool-seed42-public-verification.json); byte-identical and source-hash matched |
| Full unchanged official validator | [Passed on public revision aea6ce7](../reports/official-template-seed42-public-verification.json): fresh clone, install, two default generations and all checker gates |
| Full training-data disclosure and nonpublic-data release | [Exact prepared dataset and used nonpublic annotations public](DATA.md); owner authority confirmed; CC BY owner-controlled rights and ODbL mixed database terms retained. Historical accession/version maps and exact synthetic RL trajectories not retained; organizer interpretation not confirmed |
| Organizer acceptance of disclosed fixed-pool workflow | Not confirmed; fresh-model identity is not claimed |
| seqme / full six-category screening | Partial internal metrics only |
| Wet-lab activity, safety, synthesis/QC | Not established |
| Kaggle submission | Not performed |

Default reproduction uses a disclosed raw-candidate and score/audit cache and
reruns all selection stages. Optional `--resample` generates model candidates anew
but does not have a successful two-run certificate. Neither the public filtering
certificate nor a validator pass settles third-party data rights or organizer acceptance.
The new training release resolves the previous lack of permission/publication for
owner-controlled nonpublic data; it does not make missing provenance or activity evidence pass.

No numeric overall compliance percentage: the gates are heterogeneous and some are
blocking requirements, not interchangeable points.
