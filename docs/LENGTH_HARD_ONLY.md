# Length is a hard validity bound, not a quality target

## Policy change, 2026-09-30

The user requested reselection without penalizing peptide length beyond the
competition validity range. Length and molecular weight are now **descriptive
only**: neither enters percentile gates, the internal ranking score, or tie
breaking. Molecular weight was removed with length because it imposed a strong
indirect size filter. No length-bin quota or artificial uniformity target is added.

The sequence hard range is 8–50 aa, with standard amino acids, uniqueness,
training exact-hit exclusion and the existing complexity precheck. The existing
frozen candidate pool contains 597,555 unique sequences of 8–40 aa; reselection
does not generate or pad 41–50 aa peptides.

The five active metrics are charge at pH 7, GRAVY, aromaticity, isoelectric point
and anchor novelty. The first four retain their frozen AMP-typicality
calibration; novelty retains its high-is-better calibration. All five must pass
the same selected threshold. Physicochemical and novelty categories receive
equal ranking weight; physicochemical metrics are equally weighted within their
category. These remain a **partial-quality panel**, not the official aggregation
score or a completed six-category evaluation.

The threshold grid, CD-HIT identity 0.50 / bidirectional coverage 0.80, deterministic
cluster high/low assignment, square-root capacity allocation and tie-breaking
are unchanged. Library and Top100 retain equal high/qualified-low tiers.
Top100 reserves ten 3–4 aa motif-supported slots, five per tier, then ranks the
union by the internal score. Motifs are supporting evidence, not activity claims.

## Reselect twice from the frozen pool

Run only in an approved Lane Slurm compute allocation:

```bash
pixi run --locked select-hard-library --export \
  --work work/six-metrics-600k --output work/length-hard-only-seed42-base
pixi run --locked select-motif-top \
  --library work/length-hard-only-seed42-base --work work/six-metrics-600k \
  --output work/length-hard-only-seed42 --top-k 100 --motif-quota 10

pixi run --locked select-hard-library --export \
  --work work/six-metrics-600k --output work/length-hard-only-seed42-base-repeat
pixi run --locked select-motif-top \
  --library work/length-hard-only-seed42-base-repeat --work work/six-metrics-600k \
  --output work/length-hard-only-seed42-repeat --top-k 100 --motif-quota 10

pixi run --locked verify-length-hard-only
pixi run --locked track-length-hard-only
```

The selector defaults to seed 42, 50,000 library members and Top100. It refuses
overwriting an already completed directory. The tracking task creates one online
SwanLab audit run and sends aggregate counts and fingerprints, never protected
reference sequences. A repeated tracking invocation with the same report is not
an implicit resume.

## Candidate capacity

| Five-metric AND threshold | Candidates passing hard checks and the panel |
| --- | ---: |
| 0.50 | 48,273 |
| **0.45** | **71,592** |
| 0.40 | 95,861 |
| 0.35 | 137,543 |
| 0.30 | 166,005 |
| 0.25 | 206,089 |

The highest feasible grid value is 0.45, with high/low capacities 35,957 / 35,635.
The old seven-metric threshold of 0.30 is not used for this version.

## Verified delivery

Both independent frozen-pool selections completed and the two FASTA pairs are
byte-identical. Independent percentile/ranking recalculation and a fresh full
50,000-sequence comparison against all 40,622 known reference sequences passed.
The code suite passed 59 tests through both Pixi and uv.

| Length (aa) | Library (50,000) | Top100 |
| --- | ---: | ---: |
| 16-20 | 15,166 | 25 |
| 21-25 | 7,693 | 9 |
| 26-30 | 13,627 | 33 |
| 31-40 | 9,452 | 21 |
| 41-50 | 0 | 0 |
| 8-15 | 4,062 | 12 |

Library: mean 24.19774 aa, median 24 aa, range 8–40 aa; 27,410 VQ and
22,590 DiMA sequences, 45,568 selected clusters, 100 short-motif-supported
members. Top100: mean 24.39 aa, median 26 aa, range 9–36 aa; 46 VQ and
54 DiMA sequences, 88 clusters, 10 short-motif-supported members.

The distribution is not forced to be uniform: it reflects the existing generated
pool and the remaining quality/novelty criteria. 41–50 aa are absent from this
input pool, not rejected by a length-percentile gate.

FASTA SHA-256:

- library: `912d9a33aa2c3ff1943a52059f3bc323be3f3c9554a5399d1cd31658e5b972dc`
- top: `2354f7a6b6e3bdf4984c57c2620674cc3a7008f58d0a318386f40f0d2353c4b2`

## Verification boundary

The audit compares both complete FASTA files byte for byte, checks uniqueness,
counts, Top100 membership and tier quotas, independently recomputes the five
percentiles and category-balanced score, reconstructs Top100, checks short motif
annotations, and freshly validates the complete library against the combined
40,622-sequence official/allowed-AMP reference at Levenshtein ratio <= 0.8.

Its result is written to
[the verification report](../reports/length-hard-only-seed42-verification.json)
only on success. The two independent clustering/selection runs use the same
frozen pool. This is **not a fresh-model double-generation certificate**.

The earlier model double-run was stopped at the user's policy change; its
artifacts were retained and its GPU allocation released. The shared
`uv run generate` implementation now uses the same five quality metrics and
8–50 aa validity rule, but a full model-to-output double run under this revised
policy has not been completed. Historical seven-metric reports remain unchanged
and must not be interpreted as results for this version.

Potency/safety prediction, full embedding/collection evaluation and validated
peptide synthesis scores remain incomplete. No experimental efficacy or full
competition acceptance is claimed. This change does not publish weights, release
protected data, or submit to Kaggle.
