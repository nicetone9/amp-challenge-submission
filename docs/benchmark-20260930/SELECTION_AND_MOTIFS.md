# Selection, ranking and motif evidence

## Frozen selection

The initial existing raw pool contains 300,000 VQ and 300,000 DiMA designs; deduplication leaves 597,555. Conservative hard checks plus prechecks leave 594,243. Seven direction-normalized percentile metrics each pass >=0.30: length, charge at pH 7, GRAVY, aromaticity, isoelectric point, molecular weight (AMP typicality), and anchor novelty (higher novelty is better). Reference distributions remain fixed.

Within-category metrics receive equal weights, followed by equal weights across the available physicochemical and novelty categories. This internal ranking in [0,1] is not the official Aggregation Score. Missing potency/safety, embedding and synthesis metrics are not counted as passing.

Eligible sequences are clustered with fixed CD-HIT settings. Each cluster is split into disjoint high and qualified-low tiers; seeded singleton allocation and deterministic tie-breaking are retained. Capacity-constrained square-root cluster allocation selects 25,000 per tier for the library and 50 per tier for Top100. Top100 is drawn from the library and sorted by internal score. It is a diversity/quota-constrained ranked list, not the unconstrained global top 100 by raw score.

Library: 50,000 unique sequences, 26,861 VQ and 23,139 DiMA, spanning 38,200 selected clusters. Final motif-aware Top100: 52 VQ and 48 DiMA, spanning 81 clusters. Both selected sets have lengths 10–20 aa. This length concentration reflects the combined filters; it does not establish an optimal experimental length.

## Motif discovery and use

Exact substrings of width 3–8 were tested in eligible training AMPs. Clusters were split into 7,176 discovery and 1,809 validation sequences. Each sequence supplied one deterministic length/composition-preserving shuffle as background. Discovery tested 691,536 patterns, requiring support >=10, enrichment >=2 and BH-adjusted q<=0.05; independent-cluster validation required support >=5, enrichment >=2 and q<=0.05. Fourteen discovery hits reduced to twelve validation-supported patterns:

IWRR, IWV, IWVI, IWVIW, RIWV, RIWVI, VIW, VIWR, VIWRR, WVIW, WVIWR, WVIWRR.

All twelve are overlapping substrings of the illustrative string RIWVIWRR. They must not be interpreted as twelve independent AMP signatures. This exact-pattern set is narrow and specific to the reference data/protocol; it is not a comprehensive literature motif catalog or an independently validated AMP classifier.

The frozen motif set supports 2,183 of 597,555 unique raw candidates, 64 of 50,000 selected library members (0.128%), and 0 of the original quality-only Top100. The finalized motif-aware Top100 has 10 supported members. A match is descriptive support, not experimental evidence. Absence does not establish inactivity. The final selection reserves ten slots (five per tier) for eligible 3–4 aa motif matches, selected with the same capacity-constrained cluster allocator. The remaining ninety slots are selected from the remaining library with 45 per tier. Within each stage, the existing cluster and rank rules apply. The union is sorted by internal score. Motif support cannot rescue a failed hard constraint or percentile gate. Additional support was checked against 3,249 safe Step 1 AMP-positive training sequences intersecting the frozen reference; test, quarantine, holdout and conflict flags were excluded. This is not integration of every updated Step 1 sequence. The collective dataset (internal ID: ZYL) and OmegAMP source labels overlap and are not independent replications. The whole library additionally passes ratio <=0.8 against the 1,174 eligible AMP references absent from the official reference, covering a combined 40,622 unique sequences. No motifs were inserted into generated sequences and no candidates were manually substituted to obtain a match.

The score tables include exact match positions for all supported library members. Coordinates are zero-based, end-exclusive; repeated/overlapping matches are reported separately. This annotation does not change FASTA order or sequence content.

## Reproducibility boundary

Two independent complete selections from the same frozen pool produce byte-identical library and Top100 FASTA files, including ordering and headers. Seed 42, frozen reference/metric/tool hashes and the same execution environment are required. This does not certify fresh model regeneration, cold-clone inference, cache/no-cache or interrupted/resumed equivalence. Original checked output directories remain unchanged.
