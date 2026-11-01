#!/bin/bash
# 预算扫描各配置 + MPS chi=16/32 + PyMatching 的单进程计时（独占节点；每个 p 一个作业，约 4-5 小时）
#   TIMING_P=0.14 sbatch slurm/41_extra_timing.sbatch   -> results/paper1_extra_timing_p0.14/
#SBATCH --job-name=qecsb_extra_timing
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
source slurm/env.sh
set -euo pipefail
P=${TIMING_P:-0.10}
python scripts/run_paper1_extra.py timing --out "results/paper1_extra_timing_p${P}" --timing-p "$P" --timing-shots "${TIMING_SHOTS:-100}"
