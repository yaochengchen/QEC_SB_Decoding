#!/bin/bash
# 阶段 1 已完成、main_method.json 已写出时，从网格开始提交：
#   bash slurm/submit_from_grid.sh
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs results
test -f results/paper1_v2/main_method.json || { echo "results/paper1_v2/main_method.json 不存在"; exit 1; }
source slurm/common.sh; GRID_JOBS=${GRID_JOBS:-$N}; HARD_JOBS=${HARD_JOBS:-1}
j3=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-$((GRID_JOBS-1)) \
     --export=ALL,NSHARDS=$GRID_JOBS slurm/20_grid_run.sbatch)
j4=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --dependency=afterok:$j3 --array=0-$((HARD_JOBS-1)) \
     --export=ALL,NSHARDS=$HARD_JOBS slurm/21_grid_hardchi.sbatch)
j5=$(sbatch --parsable $EXTRA --dependency=afterok:$j4 slurm/22_grid_analyze.sbatch)
j6=$(sbatch --parsable $EXTRA --dependency=afterany:$j5 slurm/30_timing.sbatch)
echo "submitted: grid=$j3 hardchi=$j4 grid_analyze=$j5 timing=$j6"
