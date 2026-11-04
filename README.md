# SEDNA

Spherical Exploration of Diversity using Nucleic Acid sequences: every SSU rRNA sequence on one navigable globe,
from domains as continents down to species, with amplicon sequence variants (ASVs) placed on the same map.

This branch rebuilds SEDNA on [fastsne](https://github.com/rec3141/fastsne): GPU t-SNE run directly on a sphere.
The 2016 version (R `smacofSphere` MDS of the SILVA 97% seed, Cesium globe) is kept under `legacy/`.

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
| `pipeline/kmers.py` | any fasta -> k-mer table |
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
