#!/bin/bash
# 提交电路级局部性作业到 n013：  bash slurm/submit_locality.sh
#   WBIG=64 bash slurm/submit_locality.sh   （d>=11 和 bb72 的 worker 数，默认 96）
#   EXTRA="--partition=xxx" 同其他 submit 脚本
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
EXTRA=${EXTRA:-}
j=$(sbatch --parsable $EXTRA --export=ALL,WBIG=${WBIG:-96} slurm/50_circuit_locality.sbatch)
echo "submitted: locality=$j (node n013)  ->  results/circuit_locality/summary/ ; progress: results/circuit_locality/locality.log"
