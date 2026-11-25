#!/usr/bin/env bash
# A genome globe whose fragments were all embedded together (embed_all_grex.sbatch) -> viewer data at three
# densities, as the GTDB globes: an 11% and a 30% id-hash subset with per-point taxonomy (gtdb_globe_set.py +
# pack_globe.py) and every fragment with taxonomy once per genome (pack_gtdb_all.py).
#
# usage: pipeline/genome_globe_levels.sh GENOMES.tsv ALL_DIR TRAITS OUT_DIR "TITLE"
#   GENOMES.tsv  gtdb_prep.py's genome table
#   ALL_DIR      all.tax.tsv (gtdb_all_tax.py --full), all.conflicts.tsv, all.sphere.tsv, all.sculpt.tsv
#   TRAITS       traits table keyed on a GENOMES.tsv column (pack_globe.py --traits)
# writes OUT_DIR/globe{11,30} (pack_globe.py's .b64.txt files) and OUT_DIR/globe100 (.bin files)
# env: PYTHON (python3), LEVELS ("11 30": the subset percentages)
set -euo pipefail
genomes=$1; all=$2; traits=$3; out=$4; title=$5
here=$(cd "$(dirname "$0")" && pwd)
py=${PYTHON:-python3}
mkdir -p "$out"
layouts=(--layout "similarity=$all/all.sphere.tsv:$all/all.sphere.tsv" --layout "sculpt=$all/all.sculpt.tsv:$all/all.sculpt.tsv")
for f in ${LEVELS:-11 30}; do
  frac=$(awk -v f="$f" 'BEGIN {printf "%.4f", f / 100}')
  $py -I "$here/gtdb_globe_set.py" "$all/all.tax.tsv" "$genomes" "$all/all.conflicts.tsv" "$out/l$f" --extra 0 \
      --ref-shown "$frac" "${layouts[@]}"
  $py -I "$here/pack_globe.py" "$out/l$f.tax.tsv" "$out/globe$f" --layout similarity="$out/l$f.similarity.tsv" \
      --layout sculpt="$out/l$f.sculpt.tsv" --flags "$out/l$f.conflicts.tsv" --traits "$traits" --title "$title"
done
$py -I "$here/pack_gtdb_all.py" "$genomes" "$all/all.conflicts.tsv" "$traits" "$out/globe100" "${layouts[@]}" \
    --title "$title"
