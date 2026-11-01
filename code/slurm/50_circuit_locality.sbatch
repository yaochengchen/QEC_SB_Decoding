#!/bin/bash
# 第一篇附录：电路级 DEM 稳定子群的局部性（表面码 d=3..13 + bb72 Bravyi 电路）。只用 CPU，只在 n013 上跑。
#   bash slurm/submit_locality.sh            （或直接 sbatch slurm/50_circuit_locality.sbatch）
#   WBIG=64 bash slurm/submit_locality.sh     d>=11 和 bb72 的 worker 数（内存），默认 96
#SBATCH --job-name=qecsb_locality
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --nodelist=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=512
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=24:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_circuit_locality.py run --out results/circuit_locality --workers "$W" --workers-big "${WBIG:-96}"
