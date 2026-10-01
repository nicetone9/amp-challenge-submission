# Sequence-only AMP design with VQ-VAE and DiMA

Current delivery: the new seed-42 600,000-attempt pool sampled on 2026-09-30.
Older seven-metric and length-hard-only selections remain historical reports.

## Abstract

We developed two independent peptide generators using frozen residue representations
from the original ESM-2 650M model: a VQ-VAE with an autoregressive code prior and a
DiMA latent diffusion model. An equal-source pool of 600,000 attempts yielded 597,058
unique candidates. Canonical-sequence, training-overlap and conservative known-AMP
novelty checks were followed by five reference-calibrated internal metrics,
deterministic clustering and balanced high/qualified-low selection. The resulting
library contains 50,000 unique peptides and a ranked 100-member subset, with ten
short-motif quota selections. Length is only an 8–50 aa validity bound; length and
mass do not affect quality gates or ranking. The public default entry reruns
selection from this declared frozen pool. Activity, safety and synthesis success
have not been experimentally established.

## Models, data and external resources

Both routes use the original frozen `facebook/esm2_t33_650M_UR50D` representations,
without Step 0 adaptation or structural inputs. VQ-VAE uses a 512-entry,
64-dimensional codebook and a six-layer causal code prior. DiMA uses a twelve-layer
latent denoiser, 50 stochastic reverse steps and a trained ESM-initialized sequence
head. Six frozen policy checkpoints per route are cycled in a fixed order.
ANIA species-level potency and HemoPI2 hemolysis predictions were policy-reward
proxies; they are not measured MIC/HC50, MDR-isolate activity, or the final
five-metric ranking. Not every policy is claimed to improve its baseline.

Training uses the frozen September 4 milestone: 1,332,764 intact canonical sequences,
split by shared clusters into 1,066,105 train, 132,921 validation and 133,738 internal
evaluation sequences. Sampling weights are 30% AMP core, 60% AMP-source candidates
and 10% background/uncertain; background is not experimentally verified negative.
Core sources include the collective dataset used in a previous study (internal ID:
ZYL) and OmegAMP; AMPSphere provides the non-core candidate tier. Protected
test/quarantine sequences are excluded. Reference calibration uses 8,985 eligible
training AMPs; motif support uses only frozen, safe training evidence, not all later
Step 1 updates. [The exact prepared dataset and used positive annotation rows are now public](../DATA.md),
with retained processing snapshots, source fingerprints and mixed-license notices.
Owner redistribution authority is confirmed; missing historical accession/version
maps and exact original synthetic RL trajectories remain explicitly disclosed.

External resources include ESM-2, DiMA, ANIA, HemoPI2, CD-HIT and the official
antibacterial reference. ProtT5/collection and synthesizability evaluations remain
incomplete for this library. Upstream code, weights and data retain their own terms.

## Sampling, filtering and ranking

Each model contributes 300,000 raw candidates at master seed 42. Batches of 128
use architecture/batch-derived seeds; complete final batches give 600,064 physical
attempts, with deterministic tail exclusion. The training length distribution is
restricted to 8–50 aa. No peptide is truncated, manually edited or substituted.

The first full raw sampling completed, but its scoring stage failed because a
legacy descriptor guard rejected valid 41–50 aa peptides. The descriptor guard was
repaired without expanding the external-oracle domain. Raw batch counts, seeds and
hashes were verified, then CPU scoring/selection reused those completed samples
without resampling. The failed sampling run remains recorded as crashed. The
source finalization is not a successful two-run model-resampling certificate.

After deduplication, all candidates receive charge at pH 7, GRAVY, aromaticity, pI
and anchor-novelty scores. Four physicochemical descriptors use reference-calibrated
AMP typicality, while greater anchor novelty is favorable. Length and molecular
weight are descriptive only. Each of the five percentiles must pass a common
threshold; ranking weights the physicochemical mean and novelty equally.

Hard gates reject noncanonical/out-of-range sequences, training exact hits,
single-AA fraction above 0.6 and Levenshtein ratio above 0.8 against 40,622 combined
known references (39,448 official + 1,174 eligible training AMPs). Every candidate
that can pass the minimum allowed percentile is audited; the final whole library
is checked again. This is not coverage of every public AMP database.

The highest feasible threshold in the predefined 0.50/0.45/0.40/0.35/0.30/0.25
grid is 0.50: 51,041 candidates pass, with tier capacities 25,640 high / 25,401 low.
There is no fallback to a lower threshold for this delivery. CD-HIT uses identity
0.50, bidirectional coverage 0.80 and one thread with fixed ordering. Capped
square-root cluster allocation selects 25,000 per tier. The library has 46,286
selected clusters; final source proportions are not forced to stay equal.

Top100 uses the same cluster/tier rules: five 3–4 aa motif quota selections per tier,
then 45 remaining selections per tier, sorted by internal score with stable
tie-breaking. It is a diversity/quota-constrained ranking, not the unconstrained
global top 100. The library has 104 short-motif-supported members; Top100 has ten.
Frozen motif enrichment is supporting evidence, not proof of AMP function.
The internal score is not the official competition Aggregation Score.

## Outputs

| Output | Library | Ranked Top100 |
| --- | ---: | ---: |
| Unique members | 50,000 | 100, all from the library |
| VQ / DiMA | 27,360 / 22,640 | 54 / 46 |
| High / qualified-low | 25,000 / 25,000 | 50 / 50 |
| Length range | 8–50 aa | 8–36 aa |
| Mean / median length | 24.28 / 24 aa | 24.06 / 25 aa |
| Members longer than 20 aa | 60.632% | 65% |
| Short-motif-supported | 104 | 10 |

Only 490 library members (0.98%) are 41–50 aa; none are in Top100. This follows the
sampled pool and other selection criteria, not a 40-aa cutoff or length ranking.

Source-finalization SHA256:
library `bc1835e2b2a3c2e5185840b0800cbf779e679230c8751990c3a55a79b0ad1354`;
top `264c4ec0c407aace74ffb41d7911dd8241eae695d899ce3ff0b63cd2bd94ebc9`.

## Public reproduction and participation boundary

The default `uv run generate` anonymously downloads the SHA256-verified frozen
pool named in `candidate-pool.json`. Raw FASTAs, every unique candidate's score/audit
table and frozen references are declared inputs. Percentiles, threshold feasibility,
ranking, clustering, tier quotas, Top100 and final whole-library novelty validation
are rerun. No final selected library is provisioned as input. See
[installation and verification](../FROZEN_POOL_REPRODUCIBILITY.md).

[Two public cold-clone runs passed](../../reports/frozen-pool-seed42-public-verification.json)
at revision `05a99cf633084de146c3e91a897e89161c715288`: each independently downloaded
the pool, used a new locked runtime and empty uv environment/cache, and ran the
full default command. Both FASTAs are byte-identical and match source-finalization
hashes. Official sequence components and the independent five-metric check passed.
[The complete unchanged official validator also passed](../../reports/official-template-seed42-public-verification.json)
on public revision `aea6ce7ab7b50bd517b7788ed30037f83bafb914`: it cloned the
public repository anew, installed dependencies, ran the default entry twice,
checked sequence/reference gates and compared both outputs byte for byte.
Certified runtime code, locks and pool metadata are unchanged across these revisions.
These checks do not establish fresh-model resampling identity, organizer acceptance
or wet-lab activity.

The repository and optional frozen inference weights are public. Original code is
MIT, vendored DiMA is MIT, ESM-derived decoder parameters retain Meta MIT, and the
official validator is BSD-3-Clause; full notices are preserved. The exact prepared
sequence dataset and used nonpublic annotations are public. Owner-controlled rights
are CC BY 4.0; mixed database rights retain ODbL and individual third-party terms.
Historical provenance gaps and original synthetic RL trajectory retention remain
limitations. Full co-authorship eligibility is not claimed. The frozen-pool workflow is disclosed,
but organizer acceptance of this workflow has not been confirmed. No Kaggle
submission has been made. Missing potency/safety, embedding/collection and synthesis
evaluation are not counted as passing.
