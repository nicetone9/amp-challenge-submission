# Six-category calibration and mixed-pool expansion

Status: **partial implementation and partial evaluation**, 2026-09-30.
This is not the competition's unpublished Aggregation Score, not wet-lab evidence,
and not a claim that 50,000 candidates meet all six categories.

The existing `generate` entry point remains the historical single-route baseline.
It does **not** implement the new strict mixed pipeline yet. Do not label its
outputs as six-category-qualified. Existing FASTAs have not been overwritten.

## Current execution

Lane CPU job 396076: pool1 / lanec2-3-40, 4 CPU, 16 GB, 4 hours.
Lane GPU job 396156: model4 / lanec2-4-24, 1 RTX A6000, 8 CPU, 64 GB, 6 hours.
The allocated GPU was checked inside the job and was idle before sampling.
No model work runs on the login node or local Mac.

Both 50,000-sequence inputs were hash-checked. They contain 100,000 unique
sequences in total. There are 34 exact training hits; these are excluded from
qualification. The frozen milestone supplies 8,985 eligible 8–40-aa train-core
AMPs: 4,497 anchors and 4,488 calibration queries, separated by shared cluster ID.
Another 227 train-core sequences are outside the supported length/alphabet range.

This is the verified **milestone train-core subset**, not a claim to include all
later Step 1 database updates. Chemical modification metadata is absent from
the sequence-only milestone. Additional data requires an audited positive-evidence,
chemical-state and protected-split manifest before use.

## Registered metrics

See `src/amp_submission/calibration.py` and the generated private
`work/six-metrics/metric-registry.json`.

| Category | Registered metrics | Current evidence |
|---|---|---|
| Physicochemical | Length, pH7 charge, GRAVY, aromaticity, pI, molecular weight | All 8,985 references and 100,000 candidates computed |
| Published predictor potency / safety adjunct | ANIA E. coli, P. aeruginosa, S. aureus native log MIC; HemoPI2 HC50 in uM | Actual two-sequence CPU integration smoke passed; full scoring pending |
| Biological embedding | ProtT5 anchor-nearest cosine distance, set FBD and RBF MMD | Actual two-reference + two-candidate smoke passed; full scoring pending |
| Synthesizability | Validated short-peptide synthesis score | **UNAVAILABLE**; no molecular SAScore substitution |
| Novelty | 1 - nearest anchor Levenshtein ratio | All calibration queries and all candidates computed against identical disjoint anchors |
| Internal diversity | Mean pairwise normalized edit distance and uniqueness | Equal-size reference-subset gate implemented and fixture-tested; full collection evaluation pending |

Six categories are mandatory. Per-sequence metrics combine by **AND** at
percentile >= 0.50, not by mean. Missing / nonfinite / out-of-range results fail.
Set metrics are separate gates and are never copied onto individual sequences.

For increasing metrics use midpoint empirical CDF; decreasing metrics use 1-CDF.
For nonmonotone physical properties use
`T(x) = 2 * min(F(x), 1-F(x))`, then calibrate T against reference T values.
The candidate-ranking score, once all required values exist, is the equally
weighted mean of category means. It is an **internal score**, not an official score.

## Preliminary intersection: strict filtering is expensive

Only the seven currently computed sequence metrics are included in this
**diagnostic**, which cannot qualify a submission.

| Metric | Candidate percentile pass rate |
|---|---:|
| Length typicality | 33.048% |
| Charge typicality | 41.029% |
| GRAVY typicality | 48.547% |
| Aromaticity typicality | 59.177% |
| pI typicality | 52.816% |
| Molecular-weight typicality | 30.461% |
| Anchor novelty | 75.012% |
| All seven, AND | **2.443% (2,443 / 100,000)** |

Reference queries themselves pass all seven at 4.835%.
The all-registered-metric result is **incomplete**, not a biological failure rate:
required predicted scores and synthesis scores are missing.

At 2.443%, roughly 2.05 million raw candidates would be needed for 50,000
seven-metric passes, before additional metrics, duplicates and collection gates.
This is a planning estimate, not a yield guarantee.

## Reproducible Pixi tasks

Run only inside an approved Slurm allocation:

```bash
pixi install --locked
pixi run test
pixi run prepare-reference --manifest /path/manifest.csv.gz --audit /path/audit.json
pixi run score-descriptors
pixi run score-novelty
pixi run score-oracles --oracles /path/external --device cuda
pixi run score-embedding --model /path/prott5_eval --device cuda
pixi run calibrate
pixi run track-audit
```

`calibrate` writes an incomplete reference snapshot while required scores are
missing. A partial snapshot is not an official frozen reference distribution.
Completed score tables bind ordered sequences, metric implementation, dependency
lock and tool/weight fingerprints. Changed inputs do not silently resume.
Reference and candidate tables always use the same score implementation.

### Expand the raw pool

```bash
pixi run expand-pool --pool-size 300000 --max-candidates 3000000 --seed 42
```

The initial inputs remain VQ 50k + DiMA 50k. New sampling is equally split.
Default batch size is 128; seed streams derive from seed, architecture and batch
index. Batch files are immutable, content-hashed, and reused on resume. Increasing
the target retains exactly the existing batch prefix. Full-batch tail attempts
are counted against the raw budget; a target whose rounded batches exceed the
maximum is rejected. For a near-three-million full-batch target use 2,999,968.

Each branch currently finishes its requested batch range serially. This task
**creates raw batches only**; it does not yet implement automatic score/select/
refill orchestration or stop when 50,000 fully qualified sequences are found.
The first actual GPU smoke generated 128 VQ + 128 DiMA sequences successfully;
the 300k expansion is the next stage. Raw count is not qualified unique count.

## Selection and collection gates

CD-HIT 4.8.1: identity .50, two-sided coverage .80, `-l 7` (keeps 8-aa
peptides), one thread, input ordered by decreasing length then sequence hash.
Tests retain 8–10-aa sequences and reproduce clusters after input reordering.

Qualified candidates only: split each cluster into high/low halves, singleton
assignment by seed+cluster hash, capacity-capped square-root Hamilton allocation.
Both library and Top use 50% high / 50% qualified low; odd counts favor high.
Top is selected from the library while preserving cluster and tier.

Set gates use up to 2,048 candidates and 128 frozen query subsets, no replacement
inside a subset. Both sides and anchor sampling use matched sizes. PCA64 is fit
only on anchors; MMD bandwidth derives only from anchor distances. Full-library
counts and uniqueness do not use subsampling.

## Remaining release gates

- Integrate all eligible later Step 1 AMP annotations and chemical-state evidence.
- Full ANIA/HemoPI2/ProtT5 scoring and reference freeze.
- Obtain a validated short-peptide synthesis implementation or retain UNAVAILABLE.
- Connect mixed-pool ingestion, full filtering, official-reference hard novelty,
  quotas, motif annotations, and collection gates end to end.
- Parameterize the competition `generate` entry point only after this pipeline
  works; legacy generation must not be misrepresented as the new method.
- Full two-run / warm-cold cache / interrupted-uninterrupted generation checks.
- New mixed library and Top100, Feishu selection update, and official-validator run.

No relaxed thresholds, duplicate padding, automatic Kaggle submission, public
weights, data-license changes or reviewer-access changes are authorized here.

## References

- [Official competition repository](https://github.com/szczurek-lab/amp-challenge-2027)
- [seqme metrics/model API](https://seqme.readthedocs.io/en/stable/api.html)

seqme's public API and AMP tutorial are implementation references, not the final
competition scoring formula. Small-molecule synthesis scores are not validated
short-peptide synthesis scores.

## Deadline execution update

The target remains 50,000 **qualified** unique sequences, not just 50,000 raw
samples. The allocated GPU expands the balanced raw pool in stages (300k, then
600k); the physical-attempt ceiling remains 3,000,000. The CPU reference/candidate
predictor pass runs concurrently. A deadline does not relax percentile thresholds
or make UNAVAILABLE scores pass.

Expanded-pool ingestion now verifies completed-prefix metadata, hashes, branch
counts and duplicate provenance. Use a new prepared work directory when the
input-pool fingerprint changes.

Motif discovery and independent cluster-held-out validation have run on the
audited reference: 12 reliable exact patterns, with 279 hits among the original
100,000 candidates. These are support-only annotations, not measured activity.
Use `pixi run motif-check`; no protected reference sequences are uploaded.

Oracle and ProtT5 scoring support `--track` for explicit resumable SwanLab runs.
Motif and expanded-prefix integration bring the verified CPU test count to 35.
Full GPU repeat-run FASTA equality and complete six-category acceptance remain
unverified. An eight-sequence HemoPI2 profile identified pandas-heavy composition/
transition/distribution feature extraction as the CPU bottleneck; no unverified
replacement of the published predictor is deployed.
