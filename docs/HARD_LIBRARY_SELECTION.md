# Hard-compliant library selection (2026-09-30)

This is a **partial-quality, sequence-hard-compliant research library**, not a completed six-category evaluation or an official competition score.

## Policy update

The user authorized lowering the internal percentile threshold to obtain 50,000 candidates while retaining hard constraints. The strict six-category pipeline is unchanged. This separate selector never marks missing metrics as passed.

From 597,555 unique candidates, all contain standard amino acids and have lengths 8–40. After exact-reference exclusion, the conservative whole-library official-reference Levenshtein ratio <= 0.8 gate, and existing prechecks, 594,243 remain.

The fixed threshold grid is 0.50, 0.45, 0.40, 0.35, 0.30, 0.25. Select the highest feasible grid value, retaining the per-cluster high/low quotas. Seven available metrics must each pass: length typicality, charge typicality, GRAVY typicality, aromaticity typicality, pI typicality, molecular-weight typicality, and anchor novelty. These cover only the available physicochemical and novelty panel.

| Threshold | Eligible candidates |
| --- | ---: |
| 0.50 | 14,181 |
| 0.45 | 23,346 |
| 0.40 | 31,766 |
| 0.35 | 45,212 |
| **0.30** | **58,138** |
| 0.25 | 75,474 |

## Run on an allocated Lane compute node

Use the existing frozen, hash-verified pool and reference artifacts. Do not run on the login node.

```bash
pixi run --locked select-hard-library --export \
  --work work/six-metrics-600k \
  --output work/hard-library-seed42 \
  --n-sequences 50000 --top-k 100 --seed 42
```

Omit `--export` for the audit and threshold curve only. A completed output directory cannot be overwritten. The selector uses the existing deterministic CD-HIT 0.50 identity / 0.80 bidirectional coverage protocol, retains short peptides, and uses a single thread. Both the library and Top100 retain 50% high-tier and 50% qualified low-tier sequences. Category-balanced ranking uses only the seven available metrics. Motifs break ties and provide supporting evidence only.

Outputs: `library.fasta`, `top.fasta`, `library-scores.csv`, `top-scores.csv`, `audit.json`, and `complete.json`. Full sequence outputs remain on the private compute filesystem; the repository contains aggregate verification only.

## Verified result

- Library: 50,000 unique sequences, 26,861 VQ and 23,139 DiMA; 25,000 per tier; 38,200 selected clusters.
- Top100: 52 VQ and 48 DiMA, all from the library; 50 per tier; 79 clusters.
- Selected lengths: 10–20 aa. Library: 13,738 at 8–15 aa and 36,262 at 16–25 aa.
- All selected sequences pass the cached, hash-verified full official-reference ratio <= 0.8 audit.
- Frozen motif support: 64 library sequences and **0 Top100 sequences**. No motif-based efficacy claim is made.
- Two fresh selection runs from the same frozen pool produce byte-identical library and Top100 FASTA files.
- 45 unit tests pass.

See [the machine-readable verification report](../reports/hard-library-seed42-verification.json) for hashes and exact counts.

## Important limits

This threshold is an internal selection parameter, not an official competition requirement. The final 10–20 aa range reflects the combined reference-percentile filters; it is not evidence that these lengths are optimal experimentally.

Potency/MIC and safety predictors, ProtT5 and collection-level metrics, and validated short-peptide synthesizability are incomplete. These sequences are not experimentally validated AMPs. The `passed` column means only the seven-metric panel and hard sequence gates, **not all six competition categories**. Linear, unmodified free termini specify the intended chemical design, not a measured property.

The reproducibility check covers **selection from a fixed pool**, not fresh generation, cache/no-cache equivalence, interruption/resume, or cross-hardware equivalence. The competition `uv run generate` integration and full package validation are still pending; this task does not replace that entry point. No Kaggle submission or publication of weights/data has occurred.
