#!/bin/bash
# 一致性检查（独占一个节点，并行，约 20-40 分钟；可续跑）：
#   稀疏 vs 原稠密搜索内核的候选池逐位相同；快速 O(m) SOE vs 原枚举的得分和判决相同
#   结果 results_soe_fast/check_v2/SOE_FAST_CHECK.md（第一行 ALL_OK / FAIL / INCOMPLETE）
#SBATCH --job-name=qecsb_soe_check
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=512
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=06:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_soe_fast.py check --out results_soe_fast/check_v2 --workers "$W"
