#!/bin/bash
# 计时 v2（单进程，独占节点，约 2-3 小时）：13 个预算配置用稀疏内核 + 快速 SOE 重新计时，
# sbR8S100 / saC8W100 同时计原稠密搜索和原 SOE 枚举（前后对比）；另写工作量计数 work_counts.csv
#   TIMING_P=0.14 sbatch slurm/61_soe_fast_timing.sbatch   -> results_soe_fast/timing_v2_p0.14/
#SBATCH --job-name=qecsb_soe_timing
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
python scripts/run_soe_fast.py timing --out "results_soe_fast/timing_v2_p${P}" --timing-p "$P" --timing-shots "${TIMING_SHOTS:-100}"
