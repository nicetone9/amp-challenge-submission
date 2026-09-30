# Sequence-only AMP design with VQ-VAE and DiMA

## Abstract

We developed two independent peptide generators using frozen residue representations from the original ESM-2 650M model: a VQ-VAE with an autoregressive code prior and a DiMA latent diffusion model. Policies were adapted with published potency and hemolysis predictors. An equal-source pool of 600,000 designs yielded 597,555 unique candidates. Conservative sequence and reference-novelty filters were followed by training-AMP percentile calibration, deterministic clustering, and balanced high/qualified-low selection. The final library contains 50,000 unique peptides. Its ranked 100-member subset reserves ten slots for candidates with short, reference-supported 3–4-residue motifs. Both FASTA files were identical across two independent selections from the frozen pool. These are computational designs: antimicrobial activity, safety, and synthesis success have not been experimentally established.

## Method

Both routes use the original frozen `facebook/esm2_t33_650M_UR50D` representations, without Step 0 adaptation or structural inputs. VQ-VAE maps residue embeddings to a 512-entry, 64-dimensional codebook through a length-preserving CNN. A six-layer causal Transformer samples discrete codes, which the frozen codec decodes into amino acids. DiMA samples continuous latent noise with a twelve-layer denoiser and a trained sequence decoder. The challenge sampling configuration uses 50 stochastic reverse steps, not the 500 steps in the original pilot proposal.

Separate policies address broad-spectrum, Gram-positive, Gram-negative, MDR-proxy, selectivity-proxy and joint objectives. ANIA predictions for three species and HemoPI2 HC50 predictions supply reward proxies. The MDR proxy is not an MDR-isolate predictor, and the selectivity proxy is not an experimentally measured HC50/MIC ratio. The method does not claim that every policy improved over its baseline.

## Training data and external resources

Training uses the frozen September 4, 2026 milestone: 1,332,764 eligible intact canonical sequences, split by shared clusters into 1,066,105 training, 132,921 validation and 133,738 internal-evaluation sequences. The training mixture samples 30% annotated AMP core, 60% AMP-source candidates and 10% background/uncertain sequences. Core sources are eligible ZYL antimicrobial/antibacterial positives and OmegAMP curated AMPs; the candidate tier is AMPSphere excluding core. Background sequences are not treated as experimentally verified negatives. Protected test/quarantine sequences are excluded from generator training and calibration.

The current calibration uses 8,985 eligible training AMPs. Additional Step 1 motif evidence is restricted to 3,249 safe AMP-positive training sequences intersecting that frozen reference, excluding holdout, quarantine and conflict flags. Source labels may overlap. Later Step 1 updates have not been comprehensively integrated.

The official antibacterial reference, CD-HIT, ANIA and HemoPI2 are external resources. ProtT5 evaluation remains incomplete for this final library. Exact source versions, accession lists, redistribution permissions and release of nonpublic training data are not yet complete; the accompanying data statement discloses these limitations.

## Library construction and ranking

The raw pool contains 300,000 VQ and 300,000 DiMA designs. All 597,555 unique candidates have standard amino acids and lengths 8–40 aa. Filters exclude invalid sequences, duplicates, excessive single-residue content, exact known/training matches under the frozen precheck, and similarity above 0.8 to the official reference. We conservatively apply this reference threshold to the entire library, not just Top100. The finalized library also passes the same threshold against additional eligible training AMP references: the combined audited reference has 40,622 unique sequences. This is not a claim of coverage of every public AMP database.

Seven available metrics are direction-normalized against the fixed reference: length, charge, hydropathy, aromaticity, isoelectric point, molecular weight, and anchor novelty. Nonmonotonic physical descriptors use AMP typicality rather than a larger-is-better rule. Every metric must pass its threshold. To reach 50,000 candidates, the team explicitly approved a threshold grid of 0.50, 0.45, 0.40, 0.35, 0.30 and 0.25; 0.30 was the highest feasible tested threshold, with 58,138 eligible candidates. This is an internal selection criterion, not an official competition cutoff.

Candidates are clustered at 0.50 identity and 0.80 bidirectional coverage. Deterministic capacity-constrained square-root allocation produces 25,000 high-tier and 25,000 qualified-low-tier library members. Category-balanced internal scores rank the selected candidates; no official Aggregation Score is claimed.

For Top100, ten eligible short-motif candidates are selected first, five per tier. The remaining ninety slots preserve the overall 50/50 tier balance. Both stages use the same cluster allocator, exclude already selected sequences, and apply deterministic tie-breaking. The final list is sorted by internal score. No individual sequence was manually edited or replaced to insert a motif.

## Motif evidence and novelty

Motif discovery tested exact 3–8-residue substrings using cluster-separated discovery/validation AMP sets and deterministic composition-preserving shuffled backgrounds. Both stages required enrichment of at least two and BH-adjusted q<=0.05, with minimum supports of ten and five respectively. Twelve patterns passed, but they substantially overlap; they are not twelve independent AMP signatures.

The final selector counts only 3–4-residue patterns. There are 64 supported members in the library and ten in Top100. Across the 64 supported library members, the longest observed motif match is 3 aa for 54 candidates and 4 aa for ten. Their longest continuous match to any member of the combined reference is 4–6 aa and never exceeds half the candidate length. Full-sequence similarity remains <=0.8. Motif evidence supports interpretation but does not establish AMP efficacy.

## Final outputs and reproducibility

| Output | Result |
| --- | --- |
| Library | 50,000 unique peptides; 26,861 VQ and 23,139 DiMA |
| Ranked Top100 | 52 VQ and 48 DiMA; all drawn from the library |
| Selected lengths | 10–20 aa |
| Top100 tier balance | 50 high, 50 qualified low |
| Short-motif support in Top100 | 10 peptides |
| Two independent fixed-pool selections | Both FASTA files byte-identical |
| Unit tests | 49 passed |

The default seed is 42. Reference data, metric definitions, motif set, code and environment fingerprints, batch ordering and tie-breaking are recorded. Reproduction currently covers complete selection from the same frozen candidate pool. Fresh model regeneration through the portable competition entry point, cold-clone validation, and cache/resume equivalence are not yet certified. The portable `uv run generate` entry point has not yet been integrated with this final mixed selection route.

## Participation and release status

The private repository contains inference/selection code and documentation under MIT for original code. Model weights exist on Lane but have not been distributed through the repository. Organizer read access has not been granted or verified by this work. Weight delivery and read access for `RasmusML` and `szymczakpau` remain minimum-participation release gates.

Full co-authorship requirements are not claimed: public repository/weights, complete data disclosure and permissive release of nonpublic data, and full repeated-generation verification remain outstanding. Third-party data, code and weights retain their own terms. Missing potency/safety, embedding and peptide-synthesizability evaluations are explicitly marked incomplete; the artifacts are not described as passing all six quality categories. No competition submission has been made by this workflow.
