#!/usr/bin/env bash
# Build the two mOTUs globes from mOTUs-db 4.1 (Dmitrijeva, Ruscheweyh et al., NAR 2026, doi:10.1093/nar/gkag1231):
#
#   WORK/globe_motus    one point per mOTU (species-level cluster) detected in at least one sample, four layouts
#                 named as the models globe names them (<shape>[-<model>], globe.json models), so microscape.app's
#                 Taxonomy / Similarity buttons and Model menu switch between them:
#     similarity        5-mer profile of the mOTU's ten marker genes joined by N
#     sculpt            the same, with springs toward its GTDB domain, phylum, class and order centroids
#     similarity-abund  the mOTU's relative abundance across every profiled metagenome (rows of t(sample matrix)):
#                       mOTUs that occur in the same samples sit together
#     sculpt-abund      the same, with the taxonomy springs
#   WORK/globe_samples  one point per metagenome, placed by its mOTU composition, labelled by its most abundant taxon
#                 at each rank and coloured by phylum or by study
#   WORK/globe_samples_map  the metagenomes that have coordinates (with the motus-db.org table), on a map of the Earth
#   Each directory holds viewer/globe.html's .b64.txt files and microscape.app's .bin files.
#
# usage: pipeline/build_motus.sh RAW WORK
#   RAW holds, from Zenodo (https://www.motus-tool.org/database/download.html names the records):
#     db_mOTU/ (db_mOTU.tar.gz extracted: mOTUsv4.1.db.fna.gz, *.map.tsv.gz, *.gtdb.taxonomy.*.tsv.gz,
#       mOTUsv4.1.genomes.tsv.gz, mOTUsv4.1.db.blocklist.gz; the BWA index is not needed)
#     motus41.profiles_sparse.npz, motus41.profiles_{motu,sample}_index.csv.gz, mOTUs4.1.samples.tsv.gz
#   and, for the Environment labels, motusdb_org_samples.tsv: the motus-db.org Samples table
#     (python3 -I pipeline/motus_fetch_sample_table.py RAW/motusdb_org_samples.tsv); without it the globes get no
#     environment labels
# env: FASTSNE (fastsne binary, GPU), PYTHON (python3 with numpy, scipy, scikit-learn), ANCHOR (spring weights
#   domain,phylum,class,order; build_gtdb.sh's, since mOTUs are also ~97% Bacteria and SILVA's domain pull squeezes them
#   into one cap), ANCHOR_ABUND (the same for sculpt-abund), DEDUP (fastsne --dedup, 2)
# The CPU steps take a few minutes and about 3 GB; each fastsne run is under a minute on one 20 GB GPU.
set -euo pipefail
raw=$1; work=$2
here=$(cd "$(dirname "$0")" && pwd)
fastsne=${FASTSNE:-fastsne}
py=${PYTHON:-python3}
db=$raw/db_mOTU
p=$work/motus; s=$work/samples
mkdir -p "$work"
ts() { echo "[$(date +%T)] $*"; }

ts "marker genes + taxonomy, 5-mers, traits"
[ -s "$p.tax.tsv" ] || $py -I "$here/motus_prep.py" "$db" "$p" --cog all --per-motu 1 --concat --min-len 300
[ -s "$p.fsb" ] || $py -I "$here/kmers.py" "$p.fasta" "$p.fsb" --k 5
[ -s "$p.traits.tsv" ] || $py -I "$here/motus_profile_traits.py" "$raw/motus41.profiles_sparse.npz" \
    "$raw/motus41.profiles_motu_index.csv.gz" "$raw/motus41.profiles_sample_index.csv.gz" "$p" \
    --genomes "$db/mOTUsv4.1.genomes.tsv.gz"
# log10 of the number of samples each mOTU is detected in, as a colour scheme
[ -s "$p.log_prevalence.tsv" ] || awk -F'\t' 'NR == 1 {print "organism\tlog10_samples"; next} $2 > 0 {printf "%s\t%.4f\n", $1, log($2) / log(10)}' \
    "$p.prevalence.tsv" > "$p.log_prevalence.tsv"

ts "abundance features"
[ -s "$p.abund.fsb" ] && [ -s "$s.fsb" ] || \
    $py -I "$here/motus_abund_features.py" "$raw" "$p.tax.tsv" "$p.abund" "$s"
awk -F'\t' 'NR == 1 {print "organism\tstudy"; next} {print}' "$s.study.tsv" > "$s.study_color.tsv"
# anchor rows for the profiled mOTUs (motus.abund.tax.tsv is motus.tax.tsv filtered in order)
paste <(tail -n +2 "$p.tax.tsv" | cut -f1) "$p.anchor" | \
    awk -F'\t' 'NR == FNR {if (FNR > 1) k[$1] = 1; next} ($1 in k) {print $2}' "$p.abund.tax.tsv" - > "$p.abund.anchor"

env=()
if [ -s "$raw/motusdb_org_samples.tsv" ]; then
    ts "environment labels"
    $py -I "$here/motus_env_labels.py" "$raw/motusdb_org_samples.tsv" "$raw" "$s.tax.tsv" "$p.abund.tax.tsv" "$s" "$p"
    env=(1)
fi

# The 5-mer maps embed every mOTU with markers; the globe shows those that are also profiled (the abundance rows).
# --dedup gives points whose neighbourhoods overlap heavily (chains of near-identical profiles, drawn as worms) less
# repulsion, so each chain gathers into one patch
sph=(--sphere auto --h 0.28 --seed 1 --dedup "${DEDUP:-2}")
run() { local out=$1; shift; ts "fastsne -> $out"; [ -s "$out" ] || nice "$fastsne" -o "$out" "${sph[@]}" "$@"; }
run "$p.sphere.tsv" -i "$p.fsb" --transform sqrt --pca 128 --write-knn "$p.knn"
anchor=${ANCHOR:-0.0001,0.001,0.001,0.000125}
run "$p.sculpt.tsv" -i "$p.fsb" --transform sqrt --pca 128 --anchor-file "$p.anchor" --anchor "$anchor"
# SVD scores are signed, so no per-row transform
run "$p.abund.sphere.tsv" -i "$p.abund.fsb" --transform none --pca 128
# no domain spring here: co-occurrence mixes taxa, and one Bacteria centroid pulls the whole map into a cap
run "$p.abund.sculpt.tsv" -i "$p.abund.fsb" --transform none --pca 128 --anchor-file "$p.abund.anchor" \
    --anchor "${ANCHOR_ABUND:-0,0.0002,0.0002,0.000025}"
run "$s.sphere.tsv" -i "$s.fsb" --transform none --pca 128

ts "viewer data"
[ -s "$p.conflicts.tsv" ] || $py -I "$here/conflicts.py" "$p.tax.tsv" "$p.knn" 48 "$p.conflicts.tsv"
$py -I "$here/pack_globe.py" "$p.abund.tax.tsv" "$work/globe_motus" \
    --layout similarity="$p.sphere.tsv" --layout sculpt="$p.sculpt.tsv" \
    --layout similarity-abund="$p.abund.sphere.tsv" --layout sculpt-abund="$p.abund.sculpt.tsv" \
    --models "=Marker gene 5-mers;abund=Co-occurrence" \
    --layout-label "similarity=Marker gene 5-mers|Placed by marker-gene 5-mer similarity alone." \
    --layout-label "sculpt=Marker gene 5-mers + taxonomy|Marker-gene 5-mers, with each mOTU gently pulled toward its domain, phylum, class and order." \
    --layout-label "similarity-abund=Co-occurrence|Placed by relative abundance across 115,013 metagenomes: mOTUs found in the same samples sit together." \
    --layout-label "sculpt-abund=Co-occurrence + taxonomy|Co-occurrence across metagenomes, with each mOTU gently pulled toward its domain, phylum, class and order." \
    --flags "$p.conflicts.tsv" --traits "$p.traits.tsv" --color log10_samples="$p.log_prevalence.tsv" \
    ${env:+--labels "environment=Environment=$p.env.tsv" --labels-info "environment=$p.info.tsv" --labels-tip "environment=Mostly found in"} \
    --title "mOTUs-db 4.1" \
    --about "mOTUs (species-level clusters of 3.9M genomes) from mOTUs-db 4.1, GTDB R226 taxonomy, on a sphere by spherical t-SNE (fastsne): by their marker genes' 5-mer composition, or by which metagenomes they occur in."
$py -I "$here/pack_globe.py" "$s.tax.tsv" "$work/globe_samples" \
    --layout similarity="$s.sphere.tsv" \
    --layout-label "similarity=mOTU composition|Placed by mOTU composition (Hellinger distance between relative abundance profiles)." \
    --traits "$s.traits.tsv" --color study="$s.study_color.tsv" \
    ${env:+--labels "environment=Environment=$s.env.tsv" --labels-info "environment=$s.info.tsv" --default-labels environment} \
    --title "mOTUs-db 4.1 metagenomes" \
    --about "Metagenomes profiled by mOTUs-db 4.1, placed by their mOTU composition with spherical t-SNE (fastsne). Each is labelled by its most abundant taxon at every rank."
# the samples that have coordinates, on a map of the Earth
if [ -n "${env[*]}" ]; then
    $py -I "$here/pack_globe.py" "$s.geo.tax.tsv" "$work/globe_samples_map" \
        --layout similarity="$s.geo.tsv" --geo-layout similarity \
        --layout-label "similarity=Geography|Where each metagenome was sampled (samples sharing a spot are scattered around it, more widely the more there are)." \
        --labels "environment=Environment=$s.env.tsv" --labels-info "environment=$s.info.tsv" --default-labels environment \
        --traits "$s.traits.tsv" --color study="$s.study_color.tsv" \
        --title "mOTUs-db 4.1 metagenomes with coordinates" \
        --about "The metagenomes of mOTUs-db 4.1 that have coordinates (motus-db.org), on the Earth where they were sampled."
fi
for g in globe_motus globe_samples ${env:+globe_samples_map}; do
    cp "$here/../viewer/globe.html" "$work/$g/index.html"
    for f in "$work/$g"/*.b64.txt; do base64 -d "$f" > "${f%.b64.txt}.bin"; done
done
ts "done: (cd $work/globe_motus && python3 -m http.server 8765), then http://127.0.0.1:8765/"
