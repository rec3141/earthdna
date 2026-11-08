#!/bin/bash
#SBATCH --job-name=evo2_ssu
#SBATCH --account=def-rec3141
#SBATCH --partition=lgpu
#SBATCH --nodes=1
#SBATCH --gpus=l40s:2
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out

# Evo 2 7B embeddings of a FASTA on Grex: one model per L40S, the chunks split between the two GPUs.
# Resumable: resubmit the same command after a timeout and finished chunks are skipped.
# usage: sbatch embed_evo2_grex.sbatch WORK_DIR in.fasta OUT_PREFIX [embed_evo2.py options...]
#   WORK_DIR holds embed_evo2.py and the FASTA; outputs land there. After all chunks: --assemble (CPU, no GPU).
set -euo pipefail
work=$1; fasta=$2; out=$3; shift 3
type module >/dev/null 2>&1 || { set +u; source /etc/profile; set -u; }  # sbatch from a non-login shell has no module function
module load singularity
REF=$HOME/scratch/refdbs/evo2
SIF=$REF/containers/evo2.sif
export SINGULARITYENV_TRITON_CACHE_DIR=/work/.triton
export SINGULARITYENV_TRITON_LIBCUDA_PATH=/.singularity.d/libs   # NGC image's libcuda path is missing in the SIF

MODEL=${MODEL:-evo2_7b}   # evo2_20b also fits one 48 GB L40S (keep --tokens near one window); split over two GPUs it returns NaN
dd if=$REF/models/$MODEL/$MODEL.pt of=/dev/null bs=64M 2>&1 | tail -1   # warm page cache (mmap load on Lustre is slow)

n=$(grep -c '^>' "$work/$fasta")
nchunk=$(( (n + 2047) / 2048 )); half=$(( (nchunk + 1) / 2 ))
run() {  # $1 GPU, $2 first chunk, $3 end chunk, then embed_evo2.py options
  local g=$1 s=$2 e=$3; shift 3
  CUDA_VISIBLE_DEVICES=$g singularity exec --nv -B "$REF":/refdbs:ro -B "$work":/work --pwd /work "$SIF" \
    python -u embed_evo2.py "$fasta" "$out" --model_name $MODEL --local_path /refdbs/models/$MODEL/$MODEL.pt \
      --start $s --end $e "$@" > "$work/$out.gpu$g.log" 2>&1
}
a=("$@")
run 0 0 $half "${a[@]}" &
run 1 $half $nchunk "${a[@]}" &
wait
singularity exec -B "$work":/work --pwd /work "$SIF" python embed_evo2.py "$fasta" "$out" --assemble "${a[@]}"
