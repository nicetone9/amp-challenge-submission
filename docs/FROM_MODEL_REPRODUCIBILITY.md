# From-model reproducibility

The explicit fresh-model entry point is `uv run generate --resample` with `--arch mixed` and seed 42.
It samples both trained models afresh. It does **not** load an existing candidate
library, a previous raw pool, or cached candidate scores.

This integration is distinct from the frozen 2026-09-30 motif10 delivery:
that delivery demonstrated fixed-pool selection reproducibility only.
The new full-scale certificate must be produced by two successful default runs;
unit tests and small real-model smoke runs are not substitutes.

## Inputs and supported environment

Provision the private weights and their manifest under `checkpoint/`, plus the
reference-only bundle under `checkpoint/generation-reference/`.
The bundle contains AMP calibration scores, disjoint reference anchors, the
combined known-AMP novelty reference, training-sequence hashes and frozen motifs.
It contains **no generated candidates**. Public distribution of these optional model inputs is still pending.

Use Linux x86-64, Python 3.11, the checked-in Pixi and uv locks, CUDA, and CD-HIT
4.8.1 supplied by Pixi. The default mixed route does not invoke the legacy
ANIA/HemoPI2 runtime: it uses the explicitly disclosed five-metric partial panel.
The legacy single-model modes retain their external predictor prerequisites.

After frozen assets have been provisioned, install the minimal native runtime.
This creates Python, uv and CD-HIT, without R or the legacy oracle dependencies:

```bash
pixi install --locked -e runtime
pixi run -e runtime uv sync --frozen
pixi run -e runtime uv run generate --resample
```

The actual competition command is:

```bash
uv run generate --resample
```

All experiments and installation for this project must run inside approved Lane
Slurm compute nodes. No login-node or local-Mac computation is permitted.

## Frozen defaults and selection

- Seed 42; VQ and DiMA each freshly sample 300,000 candidates (600,000 total).
- Batch size 128; full final batches mean 600,064 physical attempts. The tail is
  deterministically excluded, rather than changing the last batch shape.
- Six frozen RL policy checkpoints per architecture are cycled in fixed order.
- Complete sequences are sampled from the frozen training length distribution,
  restricted to the 8–50 aa validity range. No substring truncation or sequence editing.
- All unique candidates receive the same four physicochemical quality descriptors
  and anchor-novelty score used for reference calibration. Length and molecular
  weight are recorded descriptively, not used for gates or ranking.
- Hard gates reject training exact hits, invalid sequences, excessive single-AA
  fraction and Levenshtein ratio above 0.8 against the combined known-AMP reference.
  Heavy similarity checks are performed on every candidate that could pass the
  minimum permitted percentile, and repeated on the final library.
- Select the highest feasible threshold from 0.50, 0.45, 0.40, 0.35, 0.30, 0.25.
  The threshold must pass every one of the five registered partial-panel scores.
- Preserve CD-HIT 0.50 identity / 0.80 bidirectional coverage, deterministic
  high/low tiers, square-root capacity allocation and stable tie-breaking.
- Export 50,000 unique sequences with 25,000 high / 25,000 qualified-low.
- Top100 comes from this library, with 50 high / 50 low and a minimum quota of
  ten short-motif-supported sequences (five per tier), using only 3–4 aa motifs.
- Missing quota or insufficient candidates fails; no old sequences, duplicates
  or relaxed novelty thresholds are used as filler.

The five-metric internal score is **not** the official aggregation score and is
not completion of the original six-category quality evaluation. Motifs are
supporting evidence, not proof of antimicrobial activity.

## Two independent default invocations

The previous seven-metric run in `work/from-model-repro-seed42-v2` was stopped
when the user changed the selection policy. It is not a completed certificate.
The revised default model-to-output route needs a new full-scale verification;
fixed-pool reselection is a separate, explicitly narrower guarantee.


```bash
# Replace FULL_COMMIT_SHA with the audited 40-character GitHub commit.
pixi run reproduce-from-model --phase prepare --revision FULL_COMMIT_SHA --root work/from-model-repro-length-hard-only-cold
pixi run reproduce-from-model --phase run --root work/from-model-repro-length-hard-only-cold
```

Preparation clones the public repository twice at that exact commit and checks
the inference source/locks against the audited preparing checkout. Each clone gets
its own locked Pixi runtime, new uv environment and empty uv package cache.
The preparing checkout's Python/R paths are removed from the child environment.
Each install records its interpreter prefix and confirms that rpy2 was not installed.

Only manifest-listed frozen input assets are copied into each clone; they are not
symlinked, and no generated candidates or candidate scores are provisioned.
Both candidate/output directories must be empty. Because these inputs are still
privately provisioned, this test does not certify public asset delivery.

The run phase invokes exactly `uv run generate --resample` inside each clone's minimal
Pixi runtime, sequentially on the same allocated GPU. `UV_FROZEN=1` prevents dependency changes, and
`PYTHONHASHSEED`, numerical thread settings and deterministic CUDA settings are
fixed. Each model batch has a seed derived from the master seed, architecture and
batch index. Model initialization and online tracking cannot advance the batch's
sampling RNG stream.

Our verification runs explicitly enable online SwanLab tracking via
`AMP_TRACK_ONLINE=1`; organizer inference does not require SwanLab credentials.
No offline fallback is used. Tracking never contributes to selection or FASTA IDs.

The runner refuses an already-used candidate directory and never automatically
retries a failed run. Run 2 starts only if run 1 exits successfully.

## Evidence required for completion

The runner writes `verification.json` only after:

1. Both default processes exit zero.
2. Raw VQ and DiMA sequence hashes match across runs.
3. Model, reference, source, dependency and recorded hardware identities match.
4. Both `library.fasta` files and both `top.fasta` files are byte-identical,
   including headers and ordering.
5. The retained official template's quantity, alphabet, length, unique-library,
   Top100 membership, no-reference-overlap and Top100 similarity checks pass.
6. Input asset hashes remain unchanged after generation.

Each invocation also rechecks the stricter whole-library combined-reference
similarity gate before writing outputs. A failure retains its diagnostics and
does not overwrite the previous completed delivery.

This certifies the pinned supported environment, not arbitrary hardware or
library versions. It does not grant organizer access, publish weights/data,
submit Kaggle, resolve data licensing, or establish wet-lab activity.
