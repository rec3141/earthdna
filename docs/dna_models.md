# DNA models for placing sequences: what we tested and what we found

Notes from November 2026 on the sequence representations behind microscape.app/earthdna/models: 5-mer composition,
DNABERT-2, DNABERT-S, and Evo 2 at 7B, 20B and 40B parameters, plus a 10-dimensional "phylogeny subspace" fitted on
top of Evo 2. Every number below was measured in this repo's pipeline; scripts are named so each result can be rerun.

## The short version

| | best for | weakest at | cost per 4 kb window |
|---|---|---|---|
| **5-mers** | free, no GPU, reference-free; strong within a genus (fragments of one genome cluster tightly) | which genera are related; short windows (4 kb is noisy) | ~0 (CPU, thousands/s) |
| **DNABERT-S** | a drop-in alternative to 5-mers with a GPU | no better than 5-mers at any rank | ~13 ms (GV100) |
| **DNABERT-2** | nothing here | taxonomy beyond phylum | ~13 ms (GV100) |
| **Evo 2 7B** (layer 24) | runs anywhere (bf16 on A30s); best tree-distance subspace | slower; class/order behind the 20B | ~0.7 s (A30) |
| **Evo 2 20B** (layer 18) | **best neighbours at every rank**; the models globe's default | needs FP8 (L40S/H100); deep layers collapse | ~1 s (L40S, one model per GPU) |
| **Evo 2 40B** | generation; it is the reference model | taxonomy embeddings (worse than the 20B); needs 2 x H100 | ~1 s (2 x H100) |
| **7B phylogeny subspace** | tree distance on unseen orders, in 10 numbers (400x smaller) | neighbour grouping in under-sampled classes and archaea | the 7B's, plus nothing |

## Benchmark 1: neighbours on the genus set

One 4 kb window from one species of each GTDB r226 genus (28,850 windows, `pipeline/models_set.py`). For every window,
the share of its 15 nearest neighbours (after channel standardisation for language models, PCA to 128) in the same
taxon (`pipeline/knn_eval.py`). With one window per genus, genus agreement is 0 by construction: this measures
whether related genera sit together.

| features | dims | phylum | class | order | family |
|---|---|---|---|---|---|
| 5-mers (canonical, sqrt) | 512 | 0.43 | 0.33 | 0.17 | 0.09 |
| 5-mers on the full 10 kb fragments | 512 | 0.52 | 0.43 | 0.24 | 0.14 |
| DNABERT-S | 768 | 0.42 | 0.33 | 0.16 | 0.09 |
| DNABERT-2 (earlier, GenBank-only version of the set) | 768 | 0.26 | 0.16 | 0.06 | 0.03 |
| Evo 2 7B, layer 24 | 4,096 | 0.82 | 0.70 | 0.43 | 0.22 |
| Evo 2 20B, layer 12 | 8,192 | 0.48 | 0.33 | 0.09 | 0.03 |
| Evo 2 20B, layer 16 | 8,192 | 0.80 | 0.67 | 0.39 | 0.20 |
| **Evo 2 20B, layer 18** | 8,192 | **0.87** | **0.78** | **0.53** | **0.30** |
| Evo 2 20B, layer 21 | 8,192 | 0.21 | 0.12 | 0.04 | 0.02 |
| Evo 2 40B, block 20 (residual stream) | 8,192 | 0.70 | 0.58 | 0.31 | 0.15 |
| Evo 2 40B, block 22 | 8,192 | 0.21 | 0.12 | 0.04 | 0.02 |
| Evo 2 40B, block 23 (unit-length pooled) | 8,192 | 0.16 | 0.08 | 0.03 | 0.01 |

Per class (70 classes with >= 60 genera; `per_class.py` in the analysis notes) the 20B beats the 7B in almost every
class, most in the hard ones (Binatia 0.34 -> 0.69, Kiritimatiellia 0.31 -> 0.62, Chloroflexia 0.48 -> 0.73). Some
deep-branching, MAG-only lineages stay hard for everything (UBA1384 0.19, Dojkabacteria 0.41, Aenigmatarchaeia 0.44
with the 20B).

## Benchmark 2: tree distance (the phylogeny subspace)

After Goodfire's "phylogeny manifold": a linear map to 10 dimensions trained so the angle between two species matches
their GTDB bac120 tree distance. 2,400 bacterial species x 24 windows of 4 kb; 504 species in 30 whole orders held out.
Pearson r between predicted and true tree distance on pairs touching held-out species.

| features | raw angle r | 10-d subspace, train r | **held-out r** | held-out rho |
|---|---|---|---|---|
| 5-mers | 0.10 | 0.79 | 0.62 | 0.59 |
| DNABERT-2 | 0.11 | 0.75 | 0.48 | 0.45 |
| DNABERT-S | 0.13 | 0.83 | 0.52 | 0.50 |
| Evo 2 7B layer 24 | 0.31 | 0.90 | 0.76 | 0.73 |
| Evo 2 7B, PCA-256 first | 0.23 | 0.87 | **0.76** | 0.73 |
| Evo 2 20B layer 18 | 0.31 | 0.90 | 0.64 | 0.65 |
| Evo 2 20B, PCA-256 first | 0.21 | 0.87 | 0.73 | 0.70 |

## Model by model

### 5-mers
- **Strengths.** Free, deterministic, no reference database, no GPU. Fragments of one genome cluster tightly, so on a
  globe dense with references (the GTDB globe: many fragments per genome) neighbours are right down to genus.
- **Weaknesses.** Composition converges between unrelated lineages (GC content dominates), so *between* genera it is
  weak: on one-window-per-genus sets it is far behind Evo 2. Short windows are noisy: 4 kb spreads ~4,000 counts over
  512 canonical 5-mers (~8 per k-mer); 10 kb fragments gain 0.08-0.10 at every rank.
- **Peculiarities.** Under taxonomic springs (sculpting) on a sparse set they squeeze into a band: there is little
  structure for the springs to organise.

### DNABERT-S (`zhihan1996/DNABERT-S`)
- **Strengths.** Small (117M), fast (~75 windows/s on a GV100), runs locally with `pipeline/embed_dnabert.py`. Trained
  so fragments of one species embed together.
- **Weaknesses.** Ties 5-mers at every rank; in the tree-distance test it fits training orders better than 5-mers but
  generalises worse (held-out r 0.52 vs 0.62). Its species-contrastive training rewards what composition already does.
- **Peculiarities.** Needs transformers 4.29.2; its bundled Triton flash-attention needs an old Triton and is switched
  off; mixed precision through `torch.autocast` (casting the model to fp16 breaks its ALiBi bias).

### DNABERT-2 (`zhihan1996/DNABERT-2-117M`)
- **Weaknesses.** Mean-pooled embeddings carry little taxonomy beyond phylum; worst subspace generalisation. Left off
  the globe.

### Evo 2 7B
- **Strengths.** The most tree-like representation: best held-out tree-distance correlation (0.76). Runs on GPUs without
  FP8 (A30, compute capability 8.0) through evo2's bf16 fallback (`--no-fp8`). ~1.4 windows/s per A30.
- **Weaknesses.** Behind the 20B on neighbours at class and below. Layer 28 is much worse than 24 (the useful
  information peaks mid-network); combining 24 and 28 is worse than 24 alone.
- **Peculiarities.** Second-half ("tail") pooling is slightly better than whole-window mean; reverse-strand averaging
  helps on reads. Long padded batches can fail in cuFFT with `CUFFT_INTERNAL_ERROR` instead of an out-of-memory error:
  `embed_evo2.py` halves the batch and clears the cuFFT plan cache.

### Evo 2 20B
- **Strengths.** Best neighbours at every rank (class 0.78 vs 0.70 for the 7B, order 0.53 vs 0.43).
- **Weaknesses.** Needs FP8 through Transformer Engine (compute capability >= 8.9: L40S, H100; not A30). Its tree-distance
  subspace overfits (8,192 dims -> 82k weights on 2,400 species): held-out r 0.64, or 0.73 after PCA-256, still below
  the 7B. Layers past ~20 of 24 collapse (layer 21: 0.21 at phylum).
- **Peculiarities.** The 40B minus its late layers (see the 40B section), so its layer 21 already sits in the
  high-norm region. Arc says it needs a Hopper GPU for correctness (FP8 through Transformer Engine); our 20B runs were on
  L40S (Ada), whose FP8 GEMMs mmaemu found to match Hopper's for these models, but we have not checked the 20B's
  reference loss on L40S. Run one model per GPU (fits a 48 GB L40S; ~1 window/s each), faster than splitting one model
  across two GPUs. Split across GPUs it returned NaN until vortex's opt-in Triton kernels were disabled
  (`embed_evo2.py` now enables them only on a single GPU).

### Evo 2 40B
> **Read the release notes first:** [evo2 v0.5.0](https://github.com/ArcInstitute/evo2/releases/tag/v0.5.0) explains
> that Evo 2 20B *is* the 40B with its late, high-norm layers removed ("model surgery" guided by a logit lens; no
> retraining). Everything in the Peculiarities below is that same phenomenon, rediscovered the slow way. It also means
> the 20B and 40B rows in Benchmark 1 are one model read at different points: the 20B's `mlp.l3` at layer 18 beat the
> 40B's residual stream at block 20 because of *where* it was read, not which model.

- **Strengths.** It is the model the others are measured against: on 2 x H100 our run reproduces evo2's reference test
  exactly (loss 0.2160, accuracy 91.67%). Generation through NVIDIA's free hosted API (~36 ms per nucleotide).
- **Weaknesses.** Worse than the 20B for taxonomy embeddings: its residual stream jumps to ~1e16 at block 23 (every
  position, every channel) and stays there, so later blocks add nothing measurable in bf16 and only the vector's
  direction carries information; that direction is a poor taxonomic embedding. The best usable layer (block 20, before
  the jump) gives phylum 0.70. Needs 2 x H100 80 GB (or one H200).
- **Peculiarities.**
  - Past block 22 the `mlp.l3` output (the layer we read for 7B / 20B) is exactly zero; read block outputs instead.
  - Not a loading bug: feeding the checkpoint's stored FP8 scales back in (`pipeline/evo2_fp8_scales.py`, after mmaemu)
    with an fp32 GEMM breaks the model (loss 1.40, mmaemu's figure for non-Hopper FP8). These checkpoints depend on
    Hopper's exact FP8 GEMM; on H100 the stock Transformer Engine path is the right one.
  - Store intermediate embeddings as float32: float16 overflows at 65,504.

### The 10-d phylogeny subspace
- **Strengths.** 10 numbers per window (400x less than the 7B's 4,096), and the best predictor of tree distance on
  orders it never saw.
- **Weaknesses.** It optimises distances, which are dominated by deep splits, so it blurs sister classes: on the genus
  set its class agreement is 0.43 against 0.76 for the full 7B (top 40 classes). Classes with few training species
  collapse (Campylobacteria 0.96 -> 0.11, Vulcanimicrobiia 0.98 -> 0.04), and archaea, absent from training, scatter.
- **Use it as** a compact distance estimator, not a neighbour map. A balanced training set (archaea included, several
  species per small class) is the obvious next step.

## Things that affect every Evo 2 result

- **Training-set overlap.** Evo 2 trained on GTDB r214.1 representatives plus genomes new in r220. Most GTDB r226
  representatives (our sets) were training data. Held-out tests should use genomes new in r226.
- **Taxonomy tags.** Training sequences were prefixed with `|D__...;P__...;C__...;O__...;F__...;G__...;S__...|`
  (uppercase). The model conditions on a tag; it does not write one after DNA (a window followed by `|D__` continues as
  DNA). Embeddings may partly reflect learned lineage labels rather than emergent phylogeny.
- **Genome context does not identify a genome's own 16S.** Given 20 kb of a genome new in r226, the 40B ranked that
  genome's own 16S first among 21 candidates in 4 of 46 cases (chance ~2): no MAG-to-16S link at this context length
  (`pipeline/ssu_link_set.py`, `ssu_link_score.py`).
- **Check for zeros, not only NaN.** Two failure modes here produced finite but empty vectors.

## Where each runs

| model | hardware used | throughput (4 kb windows) | job script |
|---|---|---|---|
| 5-mers | any CPU | thousands/s | `kmers.py` |
| DNABERT-S / -2 | GV100 (local) | ~75/s | `embed_dnabert.py` (uv, no install) |
| Evo 2 7B | Grex A30, bf16 | ~1.4/s per GPU | `embed_evo2_grex.sbatch --no-fp8` |
| Evo 2 20B | Grex L40S, one model per GPU | ~1.0/s per GPU | `MODEL=evo2_20b embed_evo2_grex.sbatch` |
| Evo 2 40B | Alliance H100 (Nibi, Fir), 2 per model | ~1.0/s per pair | `MODEL=evo2_40b embed_evo2_alliance.sbatch` |
| Evo 2 7B / 40B hosted | NVIDIA API (free, 40 requests/min) | minutes per embedding call | see `nvidia_evo2_api.md` |

The hosted API suits trying out generation; it is unusable for embeddings at scale (one 4 kb 7B forward pass took
265 s and returned 45 MB).
