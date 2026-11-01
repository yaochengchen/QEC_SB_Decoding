#!/bin/bash
# 通用作业数组：运行 $CMD --shard <i> --nshards $NSHARDS（由 submit_overnight.sh 设置 CMD）
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
echo "CMD: $CMD"
eval "$CMD --workers $W --shard $SLURM_ARRAY_TASK_ID --nshards $NSHARDS"
