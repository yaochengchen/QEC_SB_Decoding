#!/bin/bash
# 在网格（results/paper1_grid，含 hard/）跑完之后提交：  bash slurm/submit_extra.sh
#   解码(4 节点数组) -> 计时(独占 1 节点, 解码结束后才开始) -> 分析
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
test -d results/paper1_grid/chunks || { echo "results/paper1_grid/chunks 不存在"; exit 1; }
source slurm/common.sh; JOBS=${JOBS:-$N}
j1=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-$((JOBS-1)) --export=ALL,NSHARDS=$JOBS slurm/40_extra_run.sbatch)
j2=$(sbatch --parsable $EXTRA --dependency=afterany:$j1 slurm/41_extra_timing.sbatch)
j3=$(sbatch --parsable $EXTRA --dependency=afterok:$j1,afterany:$j2 slurm/42_extra_analyze.sbatch)
echo "submitted: extra=$j1 timing=$j2 analyze=$j3"
