#!/bin/bash
# 计时：独占一个节点（--exclusive），单进程依次测各个 d，结果写 results/paper1_timing/
#SBATCH --job-name=qecsb_timing
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_grid.py timing --out results/paper1_timing --timing-shots 100
