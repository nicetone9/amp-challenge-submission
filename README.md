# AMP Challenge — VQ-VAE and DiMA

## Frozen library and write-up (2026-09-30)

The latest private candidate delivery is the unchanged 50,000-member library plus
the short-motif-aware Top100: [write-up](docs/benchmark-20260930/WRITEUP.md),
[requirement checklist](docs/benchmark-20260930/REQUIREMENTS_CHECKLIST.md), and
[two-run verification](reports/motif-library-seed42-verification.json).

Private delivery candidate, 2026-09-29. **Not a completed competition submission.**
No Kaggle submission or reviewer access has been granted.

Two independent sequence-only generators use frozen original
`facebook/esm2_t33_650M_UR50D` residue embeddings. VQ-VAE samples a causal code
prior and decodes its codebook; DiMA samples continuous latent noise and decodes
with a trained ESM LM head. No structural input, prediction, loss or evaluation.

## Six-category screening update

[Current protocol, partial results and release gates](docs/SIX_METRIC_PROTOCOL.md).
The new mixed-pool protocol is under implementation. The legacy `generate`
command below does not yet produce six-category-qualified candidates.

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

After approved model and oracle assets have been provisioned, the competition entry
point uses the **same implementation**:

```bash
pixi run uv sync --frozen
pixi run uv run --frozen generate
```

A pre-provisioned environment with R 4.4 and the R package kaos can instead run
`uv sync --frozen && uv run --frozen generate` directly. **uv.lock does not install R.**
Official clean-clone compatibility with this non-Python dependency remains a release gate.

Default: VQ, seed42, batch128, 50,000 unique 8–40 aa sequences, 500,000 raw attempt
cap, fixed 2,048-member ranking pool, 100 ranked candidates.
`--arch dima` selects the independent diffusion route (50 stochastic reverse steps).
Output: `generate/library.fasta` and `generate/top.fasta`.
All arguments have defaults; missing or changed assets stop before changing outputs.
Repeated generation in the same output folder is supported. Fixed seeds do not
promise byte identity across different hardware, drivers or library versions.

Weights are **not yet distributed in this repository**. See [checkpoint/README.md](checkpoint/README.md).
Real-model CPU inference tests do not replace two full GPU generations.
The unchanged official checker is retained:

```bash
pixi run verify-submission https://github.com/nicetone9/amp-challenge-submission
```

It cannot yet pass the complete cold-clone path without released assets and the
external oracle environment. Do not treat component checks as organizer acceptance.

## License and release boundaries

Original packaging code: MIT. Vendored DiMA: upstream MIT, notice retained.
Official validator: BSD-3-Clause. External predictors, model weights and datasets
retain their own terms; this MIT file does not relicense them.
Protected test/quarantine data, credentials and full training sequences are not included.

## Experimental evidence

Historical source-environment exports passed official sequence component checks
and two independent Pixi/uv generations matched their hashes.
The portable package is a new integration and needs its own full regeneration.
SwanLab: [AMP_step2challenge](https://swanlab.cn/@nicetone9/AMP_step2challenge/overview).

## Candidate review

See [Feishu review documents and full candidate tables](docs/REVIEW_LINKS.md).
