#!/bin/bash
# 补充实验：相关 BP+OSD（全网格）+ CP-SB/CP-SA 预算扫描（d=9,11,13），作业数组
#SBATCH --job-name=qecsb_extra
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
python scripts/run_paper1_extra.py run --grid results/paper1_grid --out results/paper1_extra --workers "$W" \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS"
