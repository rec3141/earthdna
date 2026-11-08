#!/usr/bin/env bash
# Build a genome-fragment globe from a GTDB release, the GTDB counterpart of build.sh, as published on
# microscape.app/earthdna/gtdb:
#   every 10 kb fragment of every species representative -> canonical 5-mer table (sharded)
#   -> a reference set (REF of fragments by id hash, at least MIN_PER_GENOME from every genome)
#   -> fastsne embeds it twice: similarity, and sculpted by springs toward domain / phylum / class / order centroids
#   -> every fragment projected onto both maps
#   -> viewer data: the REF_SHOWN hash subset of the reference + EXTRA projected fragments, with genome traits
#   -> 4 kb read models (reference fragments cut to read length, fixed at their parents' positions) for reads
# usage: pipeline/build_gtdb.sh GTDB_RELEASE_DIR WORK_DIR [BIG_DIR] [RELEASE]
#   GTDB_RELEASE_DIR holds skani/genome_paths.tsv, skani/database/ (genome fna.gz), taxonomy/gtdb_taxonomy.tsv and
#   metadata/{bac120,ar53}_metadata_r<rel>.tsv.gz (https://data.gtdb.ecogenomic.org/releases/, for the traits)
#   BIG_DIR takes the all-fragment shards, models and projected maps (~150 GB for r226); default WORK_DIR
# env: FASTSNE (fastsne binary), PYTHON (python3 with numpy), REF (reference fraction, 0.30: 10.5M fragments for r226,
#   ~27 GB of GPU memory at H 0.28), MIN_PER_GENOME (5), REF_SHOWN (reference fraction drawn on the globe, 0.065),
#   EXTRA (projected fragments added to the globe, 1.5M; the viewer draws up to ~4M points), H (repulsion grid
#   spacing, 0.28; 0.38 needs about half the memory), ITERS (iterations after early exaggeration, 1750: 2000 in all),
#   ANCHOR (spring weights domain,phylum,class,order), READ_LEN (4000)
set -euo pipefail
gtdb=$1; work=$2; big=${3:-$2}; rel=${4:-226}
here=$(cd "$(dirname "$0")" && pwd)
fastsne=${FASTSNE:-fastsne}
py=${PYTHON:-python3}
mkdir -p "$work" "$big"
p="$work/gtdb$rel"; b="$big/gtdb$rel"
ts() { echo "[$(date +%T)] $*"; }

# shards of every fragment, plus the REF_SHOWN fraction as a first reference whose sequences seed the read models
ts prep
[ -s "$p.ref.fasta.gz" ] || $py "$here/gtdb_prep.py" "$gtdb" "$p" --ref "${REF_SHOWN:-0.065}" --frag 10000 --k 5 \
    --shard-dir "$big"
# the embedded reference: a bigger hash fraction (it contains the first) + a floor per genome, from the shards
ts reference
[ -s "$b.ref.fsb" ] || $py "$here/gtdb_reselect.py" "$p.all_fsb.txt" "$p.genomes.tsv" "$b" --ref "${REF:-0.30}" \
    --min-per-genome "${MIN_PER_GENOME:-5}"

common=(--transform sqrt --pca 128 --sphere 35 --h "${H:-0.28}" --iters "${ITERS:-1750}")
ts similarity
[ -s "$p.sphere.tsv" ] || "$fastsne" -i "$b.ref.fsb" -o "$p.sphere.tsv" "${common[@]}" --save-model "$b.sphere.model" \
    --write-knn "$p.knn"
# SILVA's spring weights except: domain 10x weaker (GTDB is ~97% Bacteria; 0.001 pulls nearly every point toward one
# centroid), phylum 2x and class 4x stronger (otherwise Bacillota and Clostridia break into separate patches)
ts sculpt
[ -s "$p.sculpt.tsv" ] || "$fastsne" -i "$b.ref.fsb" -o "$p.sculpt.tsv" "${common[@]}" --anchor-file "$b.ref.anchor" \
    --anchor "${ANCHOR:-0.0001,0.001,0.001,0.000125}" --save-model "$b.sculpt.model"
ts project
[ -s "$b.all.sphere.tsv" ] || "$fastsne" --project "$b.sphere.model" -i "@$p.all_fsb.txt" -o "$b.all.sphere.tsv"
[ -s "$b.all.sculpt.tsv" ] || "$fastsne" --project "$b.sculpt.model" -i "@$p.all_fsb.txt" -o "$b.all.sculpt.tsv"

ts globe
$py "$here/conflicts.py" "$b.ref.tax.tsv" "$p.knn" 48 "$p.conflicts.tsv"
$py "$here/gtdb_globe_set.py" "$b.ref.tax.tsv" "$p.genomes.tsv" "$p.conflicts.tsv" "$p.globe" --extra "${EXTRA:-1500000}" \
    --ref-shown "${REF_SHOWN:-0.065}" \
    --layout similarity="$p.sphere.tsv:$b.all.sphere.tsv" --layout sculpt="$p.sculpt.tsv:$b.all.sculpt.tsv"
$py "$here/gtdb_traits.py" "$gtdb/metadata/bac120_metadata_r$rel.tsv.gz" "$gtdb/metadata/ar53_metadata_r$rel.tsv.gz" \
    "$p.traits.tsv"
$py "$here/pack_globe.py" "$p.globe.tax.tsv" "$work/globe" --layout similarity="$p.globe.similarity.tsv" \
    --layout sculpt="$p.globe.sculpt.tsv" --flags "$p.globe.conflicts.tsv" --traits "$p.traits.tsv" \
    --title "GTDB r$rel genome fragments (10 kb)"
# the standalone viewer's blurb says "SSU rRNA sequences"; reword it in this copy only
sed -e 's/SSU rRNA sequences from/genome fragments from/' -e 's/Globe of SSU rRNA sequences/Globe of genome fragments/' \
    "$here/../viewer/globe.html" > "$work/globe/index.html"

# read models: REF_SHOWN reference fragments cut to READ_LEN (three pieces per 10 kb at 3 kb steps), each piece fixed at
# its parent's position, so reads project onto either map (kmers.py --canonical, then fastsne --project)
ts read models
r="$b.reads$((${READ_LEN:-4000} / 1000))k"
[ -s "$r.fsb" ] || $py "$here/gtdb_read_refs.py" "$p.ref.fasta.gz" "$r" --len "${READ_LEN:-4000}" --step 3000
$py "$here/gtdb_read_refs.py" "$p.ref.fasta.gz" "$r" --maps-only --layout sculpt="$p.sculpt.tsv" \
    --layout similarity="$p.sphere.tsv"
for L in sculpt similarity; do
  [ -s "$r.$L.model" ] || "$fastsne" -i "$r.fsb" -o /dev/null --transform sqrt --pca 128 --sphere 35 --fixed-y "$r.$L.tsv" \
      --save-model "$r.$L.model"
done
ts done
echo "open: (cd $work/globe && python3 -m http.server 8765) then http://127.0.0.1:8765/"
echo "place reads: kmers.py reads.fa reads.fsb --k 5 --canonical --min-len 1000; $fastsne --project $r.sculpt.model -i reads.fsb -o placed.tsv"
