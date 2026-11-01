#!/bin/bash
# d>=11 难 shot 的 MPS chi=32 复核（作业数组，量很小）
#SBATCH --job-name=qecsb_hardchi
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=512
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=06:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_grid.py hardchi --out results/paper1_grid --workers "$W" \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS"
