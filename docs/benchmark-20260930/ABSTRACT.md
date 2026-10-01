# Method abstract — current new 600k pool

We developed two independent, sequence-only peptide generators using frozen original
ESM-2 650M residue representations: a VQ-VAE with an autoregressive discrete-code
prior and a DiMA latent diffusion model. Published potency and hemolysis predictors
were policy-reward proxies, not experimental measurements. The new seed-42 pool
contains 300,000 attempts from each model and 597,058 unique candidates.
Canonical-sequence, training-overlap and conservative whole-library known-AMP
similarity checks are followed by five reference-calibrated internal metrics:
charge at pH 7, GRAVY, aromaticity, pI and anchor novelty. Every metric passes a common
percentile threshold of 0.50. Length is only an 8–50 aa validity bound; length and
mass do not contribute to quality gates or ranking. Deterministic clustering and
balanced high/qualified-low selection produce 50,000 unique peptides and a ranked
Top100 with ten short-motif quota selections. The public default entry reruns
filtering from the declared hash-verified frozen pool; model resampling is an
explicit, separately uncertified option. Full functional, embedding and synthesis
evaluation and training-data release remain incomplete. No experimental activity
or safety is claimed.

See [the current write-up](WRITEUP.md) for outputs, verification and release boundaries.
