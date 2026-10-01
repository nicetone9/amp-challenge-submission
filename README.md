# AMP Challenge — VQ-VAE and DiMA

## Current fixed-pool delivery (2026-09-30)

The current library is selected from one new seed-42 model-sampled pool:
300,000 VQ + 300,000 DiMA attempts, 597,058 unique candidates. It contains
50,000 unique peptides and a ranked Top100 drawn from that library.
The library spans 8–50 aa; Top100 spans 8–36 aa. Length and molecular weight
are descriptive only, not quality gates or ranking terms.

The default `uv run generate` downloads the hash-verified
[public pool release](https://github.com/nicetone9/amp-challenge-submission/releases/tag/frozen-pool-seed42-20260930)
named in `candidate-pool.json` and reruns filtering, clustering and selection.
It does not resample models or copy a preselected library. The candidate score/audit
table is a declared frozen input. [Two public cold-clone default runs passed](reports/frozen-pool-seed42-public-verification.json),
with byte-identical FASTAs matching the source selection.
See [reproduction](docs/FROZEN_POOL_REPRODUCIBILITY.md),
[write-up](docs/benchmark-20260930/WRITEUP.md), and
[requirements](docs/COMPETITION_REQUIREMENTS.md).

The five internal metrics are charge at pH 7, GRAVY, aromaticity, pI and anchor
novelty. The highest feasible tested common percentile threshold is 0.50.
The library has 25,000 high / 25,000 qualified-low members; Top100 has 50 per
tier, including ten short-motif quota selections. No quality or novelty gate
was relaxed to fill missing slots.

Two independent sequence-only generators use frozen original
`facebook/esm2_t33_650M_UR50D` residue embeddings. VQ-VAE samples a causal code
prior and decodes its codebook; DiMA samples continuous latent noise and decodes
with a trained sequence head. No structural input, loss or evaluation is used.
Explicit `--resample` needs CUDA and the released frozen model inputs;
[model-to-output reproducibility remains uncertified](docs/FROM_MODEL_REPRODUCIBILITY.md).

Historical September 29/30 selections and reports remain for provenance; they
must not be substituted for the current pool certificate. **No Kaggle submission
or organizer acceptance is claimed.**

## Six-category screening update

[Current protocol, partial results and release gates](docs/SIX_METRIC_PROTOCOL.md).
The default mixed generator uses the approved five-metric partial-quality policy,
not a completed six-category quality evaluation.

## Read first

- [Candidate selection guide](docs/CANDIDATE_SELECTION.md): Top100 + Top100 and both 50,000 libraries.
- [Top200 CSV](reports/top200-candidates.csv): sequences, original ranks, all proxy scores and compliance fields.
- [Exact full-library statistics](reports/library-statistics.json).
- [Competition requirement mapping](docs/COMPETITION_REQUIREMENTS.md).
- [Public training dataset, provenance and license scope](docs/DATA.md).
- [Ranking and score definitions](docs/METHOD.md).
- [Third-party notices](THIRD_PARTY_NOTICES.md).

The linked Top200 and architecture-specific libraries are historical research
artifacts, **not the current mixed-library submission files**.
No activity, MDR effect, hemolysis safety or synthesis success has been experimentally verified.
Original broad/Gram+/Gram−/MDR/selectivity/joint RL results include failures and unchanged checkpoints.

## Install and use

Linux x86-64, Python 3.11. Provision the locked native runtime with Pixi:

```bash
pixi install --locked -e runtime
pixi run -e runtime uv sync --frozen
pixi run -e runtime uv run generate
```

The default downloads and verifies the frozen candidate/reference release, then
executes filtering. CUDA, model weights and external predictor assets are not
needed for this route.

The minimal `runtime` environment supplies locked Python, uv and CD-HIT 4.8.1;
it does not contain R or the external oracle stack. Within that environment the
commands are simply `uv sync --frozen` and `uv run generate`. **uv.lock does not
install native CD-HIT.** The entry checks CD-HIT before filtering.
Legacy single-model predictor modes additionally require R, oracle assets and
`uv sync --frozen --extra legacy-oracles`; the full default Pixi environment retains them.

The frozen-pool cold-clone runner checks out a fixed GitHub commit twice, creates
a new minimal runtime and uv environment in each, and independently downloads
public hash-verified inputs. It certifies the default only after both runs match.

Default: mixed VQ + DiMA frozen pool, selection seed42, 600,000 raw candidates,
50,000 unique 8–50 aa sequences and 100 ranked candidates including ten short-motif
quota selections. All candidates receive five partial-panel scores; length and mass are descriptive only.
`--resample` explicitly enables fresh mixed sampling; `--arch vq` / `--arch dima`
also require `--resample` and retain the historical predictor-ranking routes.
Output: `generate/library.fasta` and `generate/top.fasta`.
All arguments have defaults; missing release metadata or changed assets stop before changing outputs.
Repeated generation in the same output folder is supported. Fixed seeds do not
promise byte identity across different hardware, drivers or library versions.

Frozen inference weights/reference assets are now distributed via the
[seed-42 release](https://github.com/nicetone9/amp-challenge-submission/releases/tag/frozen-pool-seed42-20260930).
See [checkpoint/README.md](checkpoint/README.md) for checksums and scope.
The default fixed-pool route needs no GPU or model weights.
The unchanged official checker is retained:

After a successful default generation, the downloaded pool contains the verified
official reference. Use a new, nonexistent checker clone directory:

```bash
pixi run -e runtime uv run python scripts/verify_submission.py \
  https://github.com/nicetone9/amp-challenge-submission \
  --dir work/official-submission-check \
  --antibacterial-fasta frozen-pool/antibacterial.fasta
```

The public cold-clone result must be checked against the current frozen-pool certificate.
The default mixed entry does not require the external oracle environment.
Do not treat component checks or private asset provisioning as organizer acceptance.

## License and release boundaries

Original packaging code: MIT. Vendored DiMA: upstream MIT, notice retained.
Official validator: BSD-3-Clause. External predictors, model weights and datasets
retain their own terms; this MIT file does not relicense them.
The exact prepared training dataset and used nonpublic positive annotations are now
[publicly released](docs/DATA.md), with [source-specific data terms](docs/DATA_LICENSING.md).
Protected test/quarantine data and credentials are not published. Historical upstream
accession/version gaps and exact synthetic RL trajectory retention remain disclosed limitations.

## Experimental evidence

Historical source-environment exports passed official sequence component checks
and two independent Pixi/uv generations matched their hashes.
The current default frozen-pool package has its own
[public cold-clone filtering certificate](reports/frozen-pool-seed42-public-verification.json).
The [complete unchanged official validator passed](reports/official-template-seed42-public-verification.json)
on public revision `aea6ce7ab7b50bd517b7788ed30037f83bafb914`: a third fresh
clone, dependency installation, two default generations, sequence/reference checks
and byte-for-byte reproducibility. Subsequent documentation-only changes leave the
certified runtime and frozen inputs unchanged.
SwanLab: [AMP_step2challenge](https://swanlab.cn/@nicetone9/AMP_step2challenge/overview).

## Candidate review

See [Feishu review documents and full candidate tables](docs/REVIEW_LINKS.md).
