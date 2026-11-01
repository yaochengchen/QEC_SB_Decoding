#!/bin/bash
# 阶段 1 解码（作业数组，每个数组任务处理 1/NSHARDS 的任务；中断后重新提交会跳过已完成的部分）
#SBATCH --job-name=qecsb_confirm
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=512
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=12:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_v2.py decode --out "${OUT:-results/paper1_v2}" --base-seed "${BASE_SEED:-2027101300}" --workers "$W" \
    --shard "$SLURM_ARRAY_TASK_ID" --nshards "$NSHARDS"
