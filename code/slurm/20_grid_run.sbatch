#!/bin/bash
# 阶段 2 网格解码（作业数组）
#SBATCH --job-name=qecsb_grid
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=512
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=24:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_grid.py run --out results/paper1_grid --workers "$W" \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS"
