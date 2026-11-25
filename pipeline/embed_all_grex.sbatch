#!/bin/bash
#SBATCH --account=def-rec3141
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=360G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
# Every fragment of a genome globe embedded at once on the CPU (fastsne_cpu --sphere, one Grex genoa node), with
# build_gtdb.sh's settings at h 0.28. GTDB r226 (35.1M fragments) takes 15-20 min per layout and ~95 GB; mOTUs-db 4.1
# (24.1M) ~10 min and ~65 GB.
#
# usage: sbatch -J NAME-sphere embed_all_grex.sbatch sphere WORK_DIR
#        sbatch -J NAME-sculpt embed_all_grex.sbatch sculpt WORK_DIR
#   WORK_DIR holds all_fsb.txt (the all-fragment shards, FSB1 or FSB2, one path per line), all.tax.tsv
#   (gtdb_all_tax.py) and all.anchor (gtdb_all_anchor.py), in the same row order. Writes all.{sphere,sculpt}.tsv and
#   .model; sphere also writes all.knn and all.conflicts.tsv.
# env: FASTSNE_CPU (fastsne_cpu; FSB2 shards need fastsne 6dfdbc5 or later), PIPELINE (this directory, for
#   conflicts.py), PY (a python3 with numpy; may be a wrapper such as a container exec), H (0.28), ITERS (1750)
set -euo pipefail
mode=$1; W=$2; cd "$W"
fs=${FASTSNE_CPU:-fastsne_cpu}
pipe=${PIPELINE:-$(dirname "$(readlink -f "$0")")}
common="--transform sqrt --pca 128 --sphere 35 --h ${H:-0.28} --iters ${ITERS:-1750} --threads ${SLURM_CPUS_PER_TASK:-$(nproc)}"
case $mode in
sphere)
  [ -s all.sphere.tsv ] || /usr/bin/time -v $fs -i @all_fsb.txt -o all.sphere.tsv $common --save-model all.sphere.model \
      --write-knn all.knn
  [ -s all.conflicts.tsv ] || /usr/bin/time -v ${PY:-python3} "$pipe/conflicts.py" all.tax.tsv all.knn 48 all.conflicts.tsv ;;
sculpt)
  [ -s all.sculpt.tsv ] || /usr/bin/time -v $fs -i @all_fsb.txt -o all.sculpt.tsv $common --anchor-file all.anchor \
      --anchor 0.0001,0.001,0.001,0.000125 --save-model all.sculpt.model ;;
*) echo "mode must be sphere or sculpt" >&2; exit 2 ;;
esac
echo done
