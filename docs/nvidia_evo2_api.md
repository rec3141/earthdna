# NVIDIA hosted Evo 2 API (build.nvidia.com), as tested 2026-11-08

Key: an `nvapi-` key from "Get API Key" on a model page, sent as `Authorization: Bearer $KEY`. Kept outside this repo
(an untracked env file); never commit it. The trial API appears to be free but rate limited (the account page says
"up to 40 rpm"); the organisation's $30/hour resource limit seems to apply to other services, not these calls.

## Endpoints

| endpoint | what | notes |
|---|---|---|
| `POST https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate` | 40B generation | body: `sequence` (uppercase DNA), `num_tokens`, `temperature` (<= 1.3), `top_k` (0-6), `top_p`, `random_seed`, `enable_logits`, `enable_sampled_probs`, `enable_elapsed_ms_per_token` |
| `POST https://health.api.nvidia.com/v1/biology/arc/evo2-7b/forward` | 7B forward pass, layer outputs | body: `sequence`, `output_layers` (1-100 names); returns `data` = base64 NPZ, `elapsed_ms` |

No hosted 40B forward (embeddings) endpoint is documented.

**Layer names:** the hosted 7B runs vortex's StripedHyena, so it wants vortex names such as `blocks.24.mlp.l3`
(what `pipeline/embed_evo2.py` uses). The documented Megatron names (`decoder.layers.24.mlp.linear_fc2`) fail with
`StripedHyena has no attribute 'decoder'`. Output: `blocks.24.mlp.l3.output`, float64, shape (1, length, 4096).

**Taxonomy prompt** (model card): optionally prepend `|D__BACTERIA;P__PSEUDOMONADOTA;C__GAMMAPROTEOBACTERIA;O__ENTEROBACTERALES;F__ENTEROBACTERIACEAE;G__ESCHERICHIA;S__ESCHERICHIA|`
(uppercase, D..S) to bias generation; it can also be placed inside the sequence. The model does not write a tag
after DNA: a window followed by `|D__` continued as DNA, and lowercase `|d__` gave a poly-A run.

## Measured

| call | server | wall | response |
|---|---|---|---|
| 40B generate 100 nt after 3 kb | 4.6 s | 5.2 s | 36 ms/token |
| 40B generate 1,000 nt after 3 kb | 37.8 s | 38.4 s | 36 ms/token |
| 7B forward 1 kb, one layer | 3.6 s | 275 s | 11 MB |
| 7B forward 4 kb, one layer | 13.3 s | 265 s | 45 MB |
| 7B forward 10 kb | - | HTTP 504 | - |

Verdict: fine for trying the 40B's generation; unusable for embeddings at scale (minutes and tens of MB per window,
against ~1 window/s per H100 with pooled 16 KB vectors on an Alliance cluster: `pipeline/embed_evo2_alliance.sbatch`).

## Training data (model card)

OpenGenome (GTDB v214.1 representatives, IMG/VR phages, IMG/PR plasmids) + 28,174 representatives new in GTDB r220 +
15,032 NCBI eukaryotic reference genomes + 41,253 metagenomes / MAGs + 32,241 organelle genomes; 9.3 T tokens.
Most GTDB r226 species representatives were therefore training data: held-out tests should use genomes new in r226.
