#!/bin/bash
# 通用单任务（分析）：运行 $CMD
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=08:00:00
source slurm/env.sh
set -euo pipefail
echo "CMD: $CMD"
eval "$CMD"
