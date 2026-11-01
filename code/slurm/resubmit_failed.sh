#!/bin/bash
# 某个数组任务超时/被杀后，只重跑这一步（已完成的块自动跳过），例如：
#   STEP=grid bash slurm/resubmit_failed.sh
set -euo pipefail
cd "$(dirname "$0")/.."
source slurm/common.sh
case "${STEP:?set STEP=confirm|grid|hardchi}" in
  confirm) N=${CONFIRM_JOBS:-1}; f=slurm/10_confirm_decode.sbatch ;;
  grid)    N=${GRID_JOBS:-$N};   f=slurm/20_grid_run.sbatch ;;
  hardchi) N=${HARD_JOBS:-1};    f=slurm/21_grid_hardchi.sbatch ;;
esac
sbatch $EXTRA --cpus-per-task=$CPUS --array=0-$((N-1)) --export=ALL,NSHARDS=$N "$f"
