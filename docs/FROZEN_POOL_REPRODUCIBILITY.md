# Frozen-pool submission reproduction

The default `uv run generate` filters a published, immutable pool of VQ and DiMA
candidates. Each branch contributes 300,000 raw sequences, for 600,000 total
attempts before deduplication. It does not sample the models again.

The release contains both raw FASTAs, the score/audit table for every unique
candidate, the frozen calibration/reference/motif inputs and SHA256 manifests.
The score table is a declared input: percentiles, quality ranking, threshold
selection, clustering, tier quotas and Top100 selection are recomputed each run.
The final whole-library known-AMP similarity check is also rerun.

## From a new clone

Use Linux x86-64 and the locked native runtime for Python 3.11, uv and CD-HIT 4.8.1:

```bash
git clone https://github.com/nicetone9/amp-challenge-submission.git
cd amp-challenge-submission
pixi install --locked -e runtime
pixi run -e runtime uv sync --frozen
pixi run -e runtime uv run generate
```

Inside that runtime the entry is simply `uv run generate`. The first invocation
downloads the release named in `candidate-pool.json`, verifies the archive and
manifest hashes, and installs it under `frozen-pool/`. Further invocations reuse
only those verified inputs and rerun filtering. Changed inputs fail.
CUDA, model weights, R and external predictor assets are not needed by this route.
uv does not install native CD-HIT; the checked-in Pixi runtime supplies it.

Outputs are `generate/library.fasta` (50,000 unique sequences) and
`generate/top.fasta` (100 members of that library). Successful repeated runs may
replace these outputs; a failed filtering run preserves the previous FASTAs.

## Unchanged filtering policy

- Canonical amino acids; length 8–50 aa as a validity bound only.
- Exclude exact training hits, excessive single-AA fraction and known-AMP
  Levenshtein ratios above 0.8.
- Five quality metrics: charge at pH 7, GRAVY, aromaticity, pI and anchor novelty.
  Length and molecular weight are descriptive and do not affect gates/ranking.
- Highest feasible common threshold from 0.50, 0.45, 0.40, 0.35, 0.30, 0.25.
- CD-HIT identity 0.50 / bidirectional coverage 0.80, one thread, stable ordering.
- Library: 25,000 high and 25,000 qualified-low, with square-root cluster quotas.
- Top100: 50 per tier, including at least five short-motif quota selections per tier.

## Verification

```bash
pixi run reproduce-frozen-pool --revision FULL_COMMIT_SHA --root work/frozen-pool-cold-reproduction
```

This creates two actual public GitHub clones at the supplied commit, each with
its own native runtime, uv environment, empty package cache and independently
downloaded public pool. It runs the default command twice, checks byte equality,
checks equality to the source pool finalization selection, and runs the retained
official validator's sequence components. A certificate is written only after
all checks pass.

This certifies filtering from the fixed pool in the pinned environment.
It does not certify repeated model sampling, complete six-category quality,
organizer acceptance or experimental activity.

## Explicit new sampling

`uv run generate --resample` retains fresh mixed-model generation for deliberate
new experiments. It requires the frozen model assets and an allocated CUDA GPU.
See [from-model reproduction](FROM_MODEL_REPRODUCIBILITY.md).
