#!/usr/bin/env bash
# Genome-fragment globe ("genome 5-mers") from the mOTUs-db 4.1 representative genomes: build_gtdb.sh's steps with
# one genome per mOTU in place of one per GTDB species.
#   every 10 kb fragment of every representative -> canonical 5-mer table (sharded)
#   -> a reference set (a hash fraction of the fragments, at least MIN_PER_GENOME from every genome)
#   -> fastsne embeds it twice: similarity, and sculpted by springs toward domain / phylum / class / order centroids
#   -> every fragment projected onto both maps
#   -> viewer data: the REF_SHOWN hash subset of the reference + EXTRA projected fragments, with the mOTU traits
#   -> 4 kb read models (reference fragments cut to read length, fixed at their parents' positions) for reads
#   The published globes instead embed every fragment together on Grex: gtdb_all_tax.py --full / gtdb_all_anchor.py on
#   WORK_DIR's all-fragment shards, embed_all_grex.sbatch, then genome_globe_levels.sh (11%, 30% and 100% globes).
#
# The genomes are presented to gtdb_prep.py as a GTDB release: WORK_DIR/release/skani/reps/<mOTU>_genomic.fna.gz
# links to REPS_DIR/<genome>.fa.gz, and release/taxonomy/gtdb_taxonomy.tsv holds each mOTU's GTDB R226 lineage from
# mOTUsv4.1.gtdb.taxonomy.80mv.tsv.gz (the labels the marker-gene globe uses). The accession is therefore the mOTU id
# (mOTUv4.0_000123, two '_' fields, as gtdb_reselect.py and gtdb_globe_set.py read it) and fragment ids are
# <mOTU>_<contig>_<start>. The representative of each mOTU is the GENOME of mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz.
#
# usage: pipeline/build_motus_genomes.sh DB_DIR REPS_DIR TRAITS WORK_DIR [BIG_DIR]
#   DB_DIR    db_mOTU/ (mOTUsv4.1.gtdb.taxonomy.{rep,80mv}.tsv.gz)
#   REPS_DIR  one <genome>.fa.gz per representative (fetch_reps.sh)
#   TRAITS    motus_profile_traits.py's traits table, keyed on the mOTU id (build_motus.sh's motus.traits.tsv)
#   BIG_DIR   takes the all-fragment shards, models and projected maps; default WORK_DIR
# env: FASTSNE, PYTHON, MIN_PER_GENOME (5), REF_SHOWN (0.065), EXTRA (1500000), ITERS (1750), ANCHOR, READ_LEN (4000)
#   REF_POINTS (6000000): the reference fraction is chosen after the fragments are counted so that the reference has
#     about this many fragments (at most 30%, as for GTDB); REF overrides it with a fixed fraction
#   FSB_U16: fastsne's tools/fsb_u16.py (default: next to FASTSNE). The k-mer tables are counts, so after the reference is
#     chosen they are rewritten as FSB2 (uint16, half the size; fastsne reads FSB2 from commit 6dfdbc5 on, not in v0.3.1) and fastsne reads
#     those; the FSB1 shards are removed once converted
#   H (0.38): repulsion grid spacing. GTDB's 10.5M-fragment reference at h 0.28 took ~27 GB of GPU memory; these
#     defaults keep the embedding inside a 20 GB card
set -euo pipefail
db=$1; reps=$2; traits=$3; work=$4; big=${5:-$4}
here=$(cd "$(dirname "$0")" && pwd)
fastsne=${FASTSNE:-fastsne}
py=${PYTHON:-python3}
mkdir -p "$work" "$big"
# p: gtdb_prep.py's outputs (its ref.* files are the REF_SHOWN subset); b: the embedded reference and its maps, under
# a different prefix so neither overwrites the other
p="$work/motusgen"; b="$big/motusgen.emb"
fsb_u16=${FSB_U16:-$(dirname "$(command -v "$fastsne")")/tools/fsb_u16.py}
ts() { echo "[$(date +%T)] $*"; }

ts release
rel="$work/release"
if [ ! -s "$rel/skani/genome_paths.tsv" ]; then
    mkdir -p "$rel/skani/reps" "$rel/taxonomy"
    zcat "$db/mOTUsv4.1.gtdb.taxonomy.80mv.tsv.gz" | awk -F'\t' 'NR > 1 {print $1 "\t" $2}' > "$rel/taxonomy/gtdb_taxonomy.tsv"
    zcat "$db/mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz" | awk -F'\t' 'NR > 1 {print $1 "\t" $2}' |
    while IFS=$'\t' read -r motu genome; do
        f="$reps/$genome.fa.gz"
        [ -s "$f" ] || { echo "missing representative genome $f" >&2; exit 1; }
        ln -sfn "$f" "$rel/skani/reps/${motu}_genomic.fna.gz"
        printf '%s\treps\n' "${motu}_genomic.fna.gz"
    done > "$rel/skani/genome_paths.tsv.part"
    mv "$rel/skani/genome_paths.tsv.part" "$rel/skani/genome_paths.tsv"
fi
echo "$(wc -l < "$rel/skani/genome_paths.tsv") genomes"

ts prep
[ -s "$p.ref.fasta.gz" ] || $py -I "$here/gtdb_prep.py" "$rel" "$p" --ref "${REF_SHOWN:-0.065}" --frag 10000 --k 5 \
    --shard-dir "$big"
total=$(awk -F'\t' 'NR > 1 {s += $5} END {print s}' "$p.genomes.tsv")
ref=${REF:-$(awk -v n="$total" -v want="${REF_POINTS:-6000000}" 'BEGIN {f = want / n; if (f > 0.30) f = 0.30; printf "%.4f", f}')}
# the shown subset is drawn from the reference, so the reference must hold it
ref=$(awk -v r="$ref" -v s="${REF_SHOWN:-0.065}" 'BEGIN {print (r < s ? s : r)}')
ts "reference: $ref of $total fragments"
[ -s "$b.ref.fsb" ] || $py -I "$here/gtdb_reselect.py" "$p.all_fsb.txt" "$p.genomes.tsv" "$b" --ref "$ref" \
    --min-per-genome "${MIN_PER_GENOME:-5}"

ts "FSB2 (uint16) copies"
[ -s "$b.ref.u16.fsb" ] || $py -I "$fsb_u16" "$b.ref.fsb" "$b.ref.u16.fsb"
if [ ! -s "$p.all_u16.txt" ]; then
    sed 's/\.fsb$/.u16.fsb/' "$p.all_fsb.txt" > "$p.all_u16.txt.part"
    paste "$p.all_fsb.txt" "$p.all_u16.txt.part" | xargs -P 6 -n 2 $py -I "$fsb_u16"
    mv "$p.all_u16.txt.part" "$p.all_u16.txt"
    xargs rm -f < "$p.all_fsb.txt"
fi

common=(--transform sqrt --pca 128 --sphere 35 --h "${H:-0.38}" --iters "${ITERS:-1750}")
ts similarity
[ -s "$p.sphere.tsv" ] || "$fastsne" -i "$b.ref.u16.fsb" -o "$p.sphere.tsv" "${common[@]}" --save-model "$b.sphere.model" \
    --write-knn "$p.knn"
ts sculpt
[ -s "$p.sculpt.tsv" ] || "$fastsne" -i "$b.ref.u16.fsb" -o "$p.sculpt.tsv" "${common[@]}" --anchor-file "$b.ref.anchor" \
    --anchor "${ANCHOR:-0.0001,0.001,0.001,0.000125}" --save-model "$b.sculpt.model"
ts project
[ -s "$b.all.sphere.tsv" ] || "$fastsne" --project "$b.sphere.model" -i "@$p.all_u16.txt" -o "$b.all.sphere.tsv"
[ -s "$b.all.sculpt.tsv" ] || "$fastsne" --project "$b.sculpt.model" -i "@$p.all_u16.txt" -o "$b.all.sculpt.tsv"

ts globe
$py -I "$here/conflicts.py" "$b.ref.tax.tsv" "$p.knn" 48 "$p.conflicts.tsv"
$py -I "$here/gtdb_globe_set.py" "$b.ref.tax.tsv" "$p.genomes.tsv" "$p.conflicts.tsv" "$p.globe" --extra "${EXTRA:-1500000}" \
    --ref-shown "${REF_SHOWN:-0.065}" \
    --layout similarity="$p.sphere.tsv:$b.all.sphere.tsv" --layout sculpt="$p.sculpt.tsv:$b.all.sculpt.tsv"
$py -I "$here/pack_globe.py" "$p.globe.tax.tsv" "$work/globe" --layout similarity="$p.globe.similarity.tsv" \
    --layout sculpt="$p.globe.sculpt.tsv" --flags "$p.globe.conflicts.tsv" --traits "$traits" \
    --title "mOTUs-db 4.1 representative genome fragments (10 kb)"
sed -e 's/SSU rRNA sequences from/genome fragments from/' -e 's/Globe of SSU rRNA sequences/Globe of genome fragments/' \
    "$here/../viewer/globe.html" > "$work/globe/index.html"

ts read models
r="$b.reads$((${READ_LEN:-4000} / 1000))k"
[ -s "$r.fsb" ] || $py -I "$here/gtdb_read_refs.py" "$p.ref.fasta.gz" "$r" --len "${READ_LEN:-4000}" --step 3000
$py -I "$here/gtdb_read_refs.py" "$p.ref.fasta.gz" "$r" --maps-only --layout sculpt="$p.sculpt.tsv" \
    --layout similarity="$p.sphere.tsv"
for L in sculpt similarity; do
  [ -s "$r.$L.model" ] || "$fastsne" -i "$r.fsb" -o /dev/null --transform sqrt --pca 128 --sphere 35 --fixed-y "$r.$L.tsv" \
      --save-model "$r.$L.model"
done
ts done
