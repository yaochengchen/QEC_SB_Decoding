#!/bin/bash
# 稀疏搜索内核 + O(m) SOE：一致性检查 + 计时 v2，两个独占节点作业并行。  bash slurm/submit_soe_fast.sh
#   结果全部写到新文件夹 results_soe_fast/（不动 results/）：
#     results_soe_fast/check_v2/SOE_FAST_CHECK.md      第一行应为 ALL_OK（约 20-40 分钟）
#     results_soe_fast/timing_v2_p0.10/TIMING_V2.md    新旧计时对比、工作量计数、CPU 型号（约 2-3 小时）
#   如果上一版（check/、timing_p0.10/）的作业还在排队或运行，可以 scancel 掉，已被 v2 取代。
#   EXTRA="--partition=xxx" / EXCLUDE= 同其他 submit 脚本；中断后重新执行即可（check 已完成的块跳过）
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs results_soe_fast
EXCLUDE=${EXCLUDE-n013}
EXTRA=${EXTRA:-}
[ -n "$EXCLUDE" ] && EXTRA="$EXTRA --exclude=$EXCLUDE"
j1=$(sbatch --parsable $EXTRA slurm/60_soe_fast_check.sbatch)
j2=$(sbatch --parsable $EXTRA --export=ALL,TIMING_P=${TIMING_P:-0.10} slurm/61_soe_fast_timing.sbatch)
echo "submitted: check=$j1  timing=$j2   ->  results_soe_fast/check_v2/ , results_soe_fast/timing_v2_p${TIMING_P:-0.10}/"
