#!/bin/bash
# 用原论文确认实验的真实 base seed 做逐 shot 复现（结果放到单独目录，不影响已有结果）：
#   BASE_SEED=<原 seed> bash slurm/submit_repro.sh
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
: "${BASE_SEED:?请设置 BASE_SEED=原确认实验的 base seed}"
export BASE_SEED OUT=${OUT:-results/paper1_v2_seed$BASE_SEED}
source slurm/common.sh
j1=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-0 --export=ALL,NSHARDS=1 slurm/10_confirm_decode.sbatch)
j2=$(sbatch --parsable $EXTRA --dependency=afterok:$j1 --export=ALL slurm/11_confirm_analyze.sbatch)
echo "submitted: decode=$j1 analyze=$j2  -> $OUT"
