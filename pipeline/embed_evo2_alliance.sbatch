#!/bin/bash
# Evo 2 embeddings on an Alliance (DRAC) H100 cluster (Nibi, Rorqual, Fir, Trillium): one model spread over the
# GPUs Slurm gives the job (evo2_40b needs 2 x H100 80 GB; evo2_7b / evo2_20b fit one). H100s run the FP8 input
# projections natively, as the 1B / 20B / 40B checkpoints were trained.
#
#   MODEL=evo2_40b sbatch --account=def-rec3141 --nodes=1 --gpus-per-node=h100:2 --cpus-per-task=12 --mem=200G \
#       --time=4:00:00 embed_evo2_alliance.sbatch WORKDIR in.fasta OUT [embed_evo2.py options, e.g. --start 0 --end 4]
#
# EVO2_ROOT holds containers/evo2.sif and models/<MODEL>/<MODEL>.pt (default /project/def-rec3141/evo2).
# Use --gpus-per-node (not --gpus): --gpus=h100:2 may give one GPU on each of two nodes. embed_evo2.py turns vortex's
# Triton kernels off on more than one GPU (they give NaN there); checked: 40B loss 0.2160 / 91.67% = evo2's reference.
#SBATCH --job-name=evo2
#SBATCH --output=%x-%j.out
set -euo pipefail
work=$1; fasta=$2; out=$3; shift 3
type module >/dev/null 2>&1 || { set +u; source /etc/profile; set -u; }
module load apptainer
ROOT=${EVO2_ROOT:-/project/def-rec3141/evo2}
MODEL=${MODEL:-evo2_40b}
export APPTAINERENV_TRITON_CACHE_DIR=/work/.triton
export APPTAINERENV_TRITON_LIBCUDA_PATH=/.singularity.d/libs   # the NGC image's libcuda path is missing in the SIF
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
dd if="$ROOT/models/$MODEL/$MODEL.pt" of=/dev/null bs=64M 2>&1 | tail -1   # warm the page cache (mmap load)
ax() { apptainer exec --nv -B "$ROOT":/refdbs:ro -B "$work":/work --pwd /work "$ROOT/containers/evo2.sif" "$@"; }
ax python -u embed_evo2.py "$fasta" "$out" --model_name "$MODEL" --local_path "/refdbs/models/$MODEL/$MODEL.pt" "$@" \
    > "$work/$out.log" 2>&1
case " $* " in *" --start "*) ;; *) ax python embed_evo2.py "$fasta" "$out" --assemble "$@" ;; esac
