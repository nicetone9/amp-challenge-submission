# AMP Challenge — VQ-VAE and DiMA

## Length-hard-only selection update (2026-09-30)

Length is now only an 8–50 aa validity bound. Length and molecular weight do not
contribute to percentile gates or ranking. See [the revised policy](docs/LENGTH_HARD_ONLY.md).
The existing frozen pool spans 8–40 aa; this update reselects from that pool and
does not invent 41–50 aa candidates. Earlier motif10 artifacts remain historical. The revised fixed-pool delivery has
passed two-run byte equality and a fresh whole-library reference audit:
[verification report](reports/length-hard-only-seed42-verification.json).

## Earlier frozen library and write-up (2026-09-30)

The earlier candidate delivery is the unchanged 50,000-member library plus
the short-motif-aware Top100: [write-up](docs/benchmark-20260930/WRITEUP.md),
[requirement checklist](docs/benchmark-20260930/REQUIREMENTS_CHECKLIST.md), and
[two-run verification](reports/motif-library-seed42-verification.json).

The default entry point now filters a frozen mixed VQ/DiMA pool with the same
partial-panel / short-motif pipeline. See [frozen-pool reproduction](docs/FROZEN_POOL_REPRODUCIBILITY.md).
New model sampling requires explicit `--resample`. The current 600,000-candidate
release and public cold-clone certificate are prepared separately from historical reports.

Candidate delivery, 2026-09-29. **Not a completed competition submission.**
No Kaggle submission has been made.

Two independent sequence-only generators use frozen original
`facebook/esm2_t33_650M_UR50D` residue embeddings. VQ-VAE samples a causal code
prior and decodes its codebook; DiMA samples continuous latent noise and decodes
with a trained ESM LM head. No structural input, prediction, loss or evaluation.

## Six-category screening update

[Current protocol, partial results and release gates](docs/SIX_METRIC_PROTOCOL.md).
The default mixed generator uses the approved five-metric partial-quality policy,
not a completed six-category quality evaluation.

## Read first

- [Candidate selection guide](docs/CANDIDATE_SELECTION.md): Top100 + Top100 and both 50,000 libraries.
- [Top200 CSV](reports/top200-candidates.csv): sequences, original ranks, all proxy scores and compliance fields.
- [Exact full-library statistics](reports/library-statistics.json).
- [Competition requirement mapping](docs/COMPETITION_REQUIREMENTS.md).
- [Training-data disclosure and unresolved permissions](docs/DATA.md).
- [Ranking and score definitions](docs/METHOD.md).
- [Third-party notices](THIRD_PARTY_NOTICES.md).

The Top200 is a research decision sheet, **not a single official Top200 submission**.
Each architecture retains its own library and ranked Top100.
No activity, MDR effect, hemolysis safety or synthesis success has been experimentally verified.
Original broad/Gram+/Gram−/MDR/selectivity/joint RL results include failures and unchanged checkpoints.

## Install and use

Linux x86-64, Python 3.11. Experiment setup and CPU tests use Pixi:

```bash
pixi install --locked
pixi run test
pixi run check-assets
pixi run generate
```

The default entry downloads and verifies the frozen candidate/reference release
named in `candidate-pool.json`, then executes the fixed filtering implementation.
Model weights and external predictor assets are only needed for explicit resampling:

```bash
pixi install --locked -e runtime
pixi run -e runtime uv sync --frozen
pixi run -e runtime uv run generate
```

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

Weights are **not yet distributed in this repository**. See [checkpoint/README.md](checkpoint/README.md).
The default fixed-pool route needs no GPU or model weights.
The unchanged official checker is retained:

```bash
pixi run verify-submission https://github.com/nicetone9/amp-challenge-submission
```

The public cold-clone result must be checked against the current frozen-pool certificate.
The default mixed entry does not require the external oracle environment.
Do not treat component checks or private asset provisioning as organizer acceptance.

## License and release boundaries

Original packaging code: MIT. Vendored DiMA: upstream MIT, notice retained.
Official validator: BSD-3-Clause. External predictors, model weights and datasets
retain their own terms; this MIT file does not relicense them.
Protected test/quarantine data, credentials and full training sequences are not included.

## Experimental evidence

Historical source-environment exports passed official sequence component checks
and two independent Pixi/uv generations matched their hashes.
The default frozen-pool package needs its own public cold-clone filtering certificate.
SwanLab: [AMP_step2challenge](https://swanlab.cn/@nicetone9/AMP_step2challenge/overview).

## Candidate review

See [Feishu review documents and full candidate tables](docs/REVIEW_LINKS.md).
