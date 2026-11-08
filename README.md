# EarthDNA

Every SSU rRNA sequence on one navigable globe,
from domains as continents down to species, with amplicon sequence variants (ASVs) placed on the same map.

EarthDNA is built on [fastsne](https://github.com/rec3141/fastsne): GPU t-SNE run directly on a sphere. It grew out of
SEDNA (Spherical Exploration of Diversity using Nucleic Acid sequences, 2016: R `smacofSphere` MDS of the SILVA 97%
seed on a Cesium globe), which is kept under `legacy/`.

## How the globe is made

1. **Features.** Each reference becomes a 5-mer count profile (1,024 values), square-root transformed and reduced
   to 128 principal components. On SILVA 144 NR99 the 15 nearest neighbours in this space share the query's genus
   87% of the time with exact search, the same as vsearch identity neighbours; longer k-mers add nothing once PCA
   is applied, and 128 components recover most of what 64 lose.
2. **Layout.** `fastsne --sphere 35` embeds all 905,628 NR99 sequences on a sphere (about 50 s on a Quadro GV100):
   t-SNE with the Cauchy kernel on chord distance, repulsion from a 3D FFT over the sphere's bounding cube.
3. **Sculpting (optional).** `--anchor` adds weak springs toward each sequence's domain, group, class and order
   centroids, which gathers domains into continents with ocean between them and phyla into countries, while
   sequence similarity still decides the detail.
4. **Label checks.** Sequences whose neighbours mostly carry another label are flagged. In SILVA 144 this finds
   about 5,600 land-plant nuclear 18S genes from genome assemblies labelled as mitochondria or chloroplasts.

## Build

```
# SILVA_144_SSURef_NR99_tax_silva_trunc.fasta.gz and tax_slv_ssu_144.txt.gz in $SILVA
FASTSNE=/path/to/fastsne pipeline/build.sh $SILVA work 144
(cd work/globe && python3 -m http.server 8765)      # then http://127.0.0.1:8765/
```

| script | does |
|---|---|
| `pipeline/silva_prep.py` | SILVA fasta + rank file -> k-mer table (`.fsb`), taxonomy table, anchor groups |
| `pipeline/conflicts.py` | kNN graph + taxonomy -> sequences whose label disagrees with their neighbours |
| `pipeline/pack_globe.py` | layouts + taxonomy -> the viewer's data files |
| `pipeline/region_refs.py` | aligned references -> one amplicon region per reference, for ASV placement |
| `pipeline/kmers.py` | any fasta -> k-mer table (`--rc` both strands, `--canonical` 512 columns as for GTDB) |
| `pipeline/build_gtdb.sh` | the GTDB genome-fragment globe and its read models, end to end |
| `pipeline/gtdb_prep.py` | GTDB genomes -> every 10 kb fragment's canonical 5-mers (sharded), reference sequences |
| `pipeline/gtdb_reselect.py` | a bigger reference from the shards: id-hash fraction + a floor per genome |
| `pipeline/gtdb_globe_set.py` | reference subset + projected fragments -> one point set for `pack_globe.py` |
| `pipeline/gtdb_traits.py` | GTDB metadata -> per-genome traits (cultured, named, type material) for `pack_globe.py --traits` |
| `pipeline/gtdb_read_refs.py` | reference fragments cut to read length, fixed at their parents' positions (read models) |
| `pipeline/render_preview.py` | a still image of a packed globe |
| `pipeline/pack_overlay.py` | placed reads (fastsne `--project`, per layout) -> an overlay the viewer draws on top |
| `pipeline/embed_evo2.py` | any fasta -> Evo 2 layer embeddings as `.fsb` (with `embed_evo2*_grex.sbatch`) |
| `pipeline/knn_eval.py` | feature tables -> how often nearest neighbours share each rank's label |
| `pipeline/pilot_sample.py`, `pilot_reads.py` | genus-stratified SSU / labelled nanopore-read sets for comparing features |
| `viewer/globe.html` | WebGL2 globe: equal-area view, layouts, labels by zoom, search, label conflicts |

Taxonomic ranks come from SILVA's rank file rather than positions in the path, because eukaryote paths stop at
whatever node SILVA reached. Species names come from the source-organism field and are kept only when their genus
agrees with SILVA's (for bacteria the organism is often the host).

## Placing ASVs

A short amplicon's k-mer profile is not comparable with a full-length one, so ASVs are placed with a model built
from the references cut to the same region, carrying the full-length map positions:

```
python3 pipeline/region_refs.py SILVA_144_SSURef_NR99_tax_silva_full_align_trunc.fasta.gz \
    GTGYCAGCMGCCGCGGTAA GGACTACNVGGGTWTCTAAT v4_refs.fasta          # 515F-Y / 806R-B
python3 pipeline/kmers.py v4_refs.fasta v4_refs.fsb
fastsne -i v4_refs.fsb -o /dev/null --transform sqrt --pca 128 --sphere 35 \
    --fixed-y work/nr99.sphere.tsv --save-model v4.model
python3 pipeline/kmers.py asvs.fasta asvs.fsb
fastsne --project v4.model -i asvs.fsb -o asvs.placed.tsv          # id lon lat nn_dist
```

## Genome globe (GTDB)

`pipeline/build_gtdb.sh` builds the same kind of globe from GTDB genomes (r226: 143,614 species representatives,
35.1M non-overlapping 10 kb fragments, canonical 5-mers). fastsne embeds a 30% reference (10.5M fragments, chosen by
id hash with at least 5 from every genome; 2000 iterations, ~27 GB of GPU memory at h 0.28) and every fragment is
projected onto it. The sculpted map uses SILVA's springs with the domain pull 10x weaker (GTDB is ~97% Bacteria),
phylum 2x and class 4x stronger, which keeps Bacillota and Clostridia whole. The viewer draws the 6.5% hash subset
of the reference plus 1.5M projected fragments (3.9M points), with GTDB traits to highlight cultured isolates, named
species and type material.

Of a reference fragment's 15 nearest 5-mer neighbours, 34% come from the same genome, 44% from other species of its
genus and 22% from other genera: species mix within a genus, genera mostly do not.

Nanopore reads are placed with 4 kb read models (reference fragments cut to read length at the parents' positions,
the genome counterpart of the amplicon region models). On three archived barcodes (3.6M reads), the 15 reference
fragments nearest a read share its sendsketch label (ANI >= 95) at phylum 0.89-0.93 and genus 0.10-0.58, depending on
how well the sample's taxa are represented in GTDB.

## Evo 2 (experimental)

`pipeline/embed_evo2.py` writes mean or second-half averaged Evo 2 activations of chosen layers (optionally averaged
with the reverse strand) as `.fsb` tables, batching by length; `embed_evo2_grex.sbatch` (7B, one model per GPU) and
`embed_evo2_20b_grex.sbatch` (20B over two L40S) run it on Grex. `knn_eval.py` compares feature tables on
genus-stratified sets. On 8,570 labelled 3-5 kb nanopore reads, 7B layer 24 (second half, both strands, channels
z-scored) put a read's 15 nearest neighbours in its order 66% of the time against 42% for canonical 5-mers, and in its
genus 21% against 15%; layer 28 was far worse than 5-mers.
