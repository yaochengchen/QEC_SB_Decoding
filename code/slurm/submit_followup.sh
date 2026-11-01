#!/bin/bash
# v4 补跑（网格 + 补充实验 + 过夜延伸都已完成后提交）：  bash slurm/submit_followup.sh
#   1) 补充实验数组（所有可用节点）：预算扫描延到 p=.18（5 个配置）+ d=11/13 的预算扩展行 CP-SB 8x400（MAP/r1/r2）
#      已有的块全部跳过；低 p 延伸格（p<=.06）默认只取前 10 万 shot（QECSB_FULL_LOWP_CAP 可改）
#   2) 计时：p=.10 与 p=.14 各一个独占节点作业（解码数组结束后开始；两个作业并行）
#   3) 补充分析（含新表 budget_scaled_*、budget_paired_matched、ler_vs_d_budget_scaled 和三张新图）
# 先看工作量：  QECSB_LOWP_EXT=1 python scripts/run_paper1_extra.py estimate --out results/paper1_extra
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
test -d results/paper1_grid/chunks || { echo "results/paper1_grid/chunks 不存在"; exit 1; }
source slurm/common.sh
export QECSB_LOWP_EXT=1                      # 让任务规划看到 p=.04 和 d>=11 的 p=.05/.06 格
X="python scripts/run_paper1_extra.py"
A="sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-$((N-1)) --export=ALL,NSHARDS=$N"
j1=$(CMD="$X run --grid results/paper1_grid --out results/paper1_extra" $A -J qecsb_followup slurm/generic_array.sbatch)
t1=$(sbatch --parsable $EXTRA --export=ALL,TIMING_P=0.10 --dependency=afterany:$j1 slurm/41_extra_timing.sbatch)
t2=$(sbatch --parsable $EXTRA --export=ALL,TIMING_P=0.14 --dependency=afterany:$j1 slurm/41_extra_timing.sbatch)
j3=$(CMD="$X analyze --grid results/paper1_grid --out results/paper1_extra" \
     sbatch --parsable $EXTRA --export=ALL -J qecsb_followup_analyze --dependency=afterok:$j1,afterany:$t1:$t2 slurm/generic_single.sbatch)
echo "submitted: extra=$j1 (array 0-$((N-1)))  timing_p0.10=$t1  timing_p0.14=$t2  analyze=$j3"
echo "中断/超时后重新执行本脚本即可续跑（已完成的块跳过）；只重做分析： CMD=\"$X analyze --grid results/paper1_grid --out results/paper1_extra\" sbatch $EXTRA --export=ALL -J qecsb_followup_analyze slurm/generic_single.sbatch"
