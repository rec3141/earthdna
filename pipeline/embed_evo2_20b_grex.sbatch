#!/bin/bash
#SBATCH --job-name=evo2_20b
#SBATCH --account=def-rec3141
#SBATCH --partition=lgpu
#SBATCH --nodes=1
#SBATCH --gpus=l40s:2
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

# Evo 2 20B embeddings of a FASTA on Grex: one model split over both L40S (it needs ~40 GB of weights).
# Resumable like embed_evo2_grex.sbatch; joins the parts when all chunks are done.
# usage: sbatch embed_evo2_20b_grex.sbatch WORK_DIR in.fasta OUT_PREFIX [embed_evo2.py options...]
#   default layers: blocks 12, 16, 18, 21 of 24 (0.5, 0.67, 0.75, 0.875 of depth); override with --layers
# WARNING: the 20B split over two GPUs returns NaN for most sequences (seen on Grex L40S, 2026-11-08). Run it one per
# GPU instead: MODEL=evo2_20b sbatch -p lgpu --gpus=l40s:2 embed_evo2_grex.sbatch ... --tokens 4096
set -euo pipefail
work=$1; fasta=$2; out=$3; shift 3
type module >/dev/null 2>&1 || { set +u; source /etc/profile; set -u; }  # sbatch from a non-login shell
module load singularity
REF=$HOME/scratch/refdbs/evo2
SIF=$REF/containers/evo2.sif
export SINGULARITYENV_TRITON_CACHE_DIR=/work/.triton
export SINGULARITYENV_TRITON_LIBCUDA_PATH=/.singularity.d/libs   # NGC image's libcuda path is missing in the SIF

dd if=$REF/models/evo2_20b/evo2_20b.pt of=/dev/null bs=64M 2>&1 | tail -1   # warm page cache (mmap load on Lustre)

sx() { singularity exec --nv -B "$REF":/refdbs:ro -B "$work":/work --pwd /work "$SIF" "$@"; }
layers=--layers=blocks.12.mlp.l3,blocks.16.mlp.l3,blocks.18.mlp.l3,blocks.21.mlp.l3
sx python -u embed_evo2.py "$fasta" "$out" --model_name evo2_20b --local_path /refdbs/models/evo2_20b/evo2_20b.pt \
    --tokens 16384 "$layers" "$@" > "$work/$out.log" 2>&1
sx python embed_evo2.py "$fasta" "$out" --assemble "$layers" "$@"
