#!/usr/bin/env bash
# Build the SSU globe from a SILVA release: k-mer table -> fastsne spherical maps -> viewer data.
# usage: pipeline/build.sh SILVA_DIR WORK_DIR [RELEASE]
#   SILVA_DIR holds SILVA_<rel>_SSURef_NR99_tax_silva_trunc.fasta.gz and tax_slv_ssu_<rel>.txt.gz
#   (https://ftp.arb-silva.de/current/Exports/ and Exports/taxonomy/)
# needs fastsne (GPU build, --sphere) on PATH or in $FASTSNE, and python3 with numpy
set -euo pipefail
silva=$1; work=$2; rel=${3:-144}
here=$(cd "$(dirname "$0")" && pwd)
fastsne=${FASTSNE:-fastsne}
mkdir -p "$work"
p="$work/nr99"

[ -s "$p.fsb" ] || python3 "$here/silva_prep.py" "$silva/SILVA_${rel}_SSURef_NR99_tax_silva_trunc.fasta.gz" \
    "$silva/tax_slv_ssu_${rel}.txt.gz" "$p" --k 5

# 5-mer counts, sqrt transform, 128 PCs; R 35 with h 0.28 keeps the repulsion error near 1% on a 512^3 grid
common=(--transform sqrt --pca 128 --sphere 35 --h 0.28)
"$fastsne" -i "$p.fsb" -o "$p.sphere.tsv" "${common[@]}" --save-model "$p.sphere.model" --write-knn "$p.knn"
# the same, with weak springs toward each domain / group / class / order centroid: continents and countries
"$fastsne" -i "$p.fsb" -o "$p.sculpt.tsv" "${common[@]}" --anchor-file "$p.anchor" --anchor 0.001,0.0005,0.00025,0.000125

python3 "$here/conflicts.py" "$p.tax.tsv" "$p.knn" 48 "$p.conflicts.tsv"
python3 "$here/pack_globe.py" "$p.tax.tsv" "$work/globe" --layout similarity="$p.sphere.tsv" \
    --layout sculpt="$p.sculpt.tsv" --flags "$p.conflicts.tsv" --title "SILVA $rel SSU Ref NR99"
cp "$here/../viewer/globe.html" "$work/globe/index.html"
echo "open: (cd $work/globe && python3 -m http.server 8765) then http://127.0.0.1:8765/"
