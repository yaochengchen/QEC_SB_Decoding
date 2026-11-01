#!/bin/bash
# 冒烟测试通过后在登录节点执行：  bash slurm/submit_all.sh
# 一次提交全部作业，用依赖关系串起来：
#   阶段1解码(数组) -> 阶段1分析+复现闸门 -> 网格(数组) -> 难shot复核(数组) -> 网格分析
#   计时作业放在最后（网格跑完、节点空出来以后）
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs results

source slurm/common.sh            # N = 可用节点数（排除 n013），CPUS、EXTRA、EXCLUDE 见 common.sh
GRID_JOBS=${GRID_JOBS:-$N}       # 网格数组大小 = 节点数
CONFIRM_JOBS=${CONFIRM_JOBS:-1}  # 阶段1只有 ~16 核时，1 个节点足够
HARD_JOBS=${HARD_JOBS:-1}

j1=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-$((CONFIRM_JOBS-1)) \
     --export=ALL,NSHARDS=$CONFIRM_JOBS slurm/10_confirm_decode.sbatch)
j2=$(sbatch --parsable $EXTRA --dependency=afterok:$j1 slurm/11_confirm_analyze.sbatch)
j3=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --dependency=afterok:$j2 --array=0-$((GRID_JOBS-1)) \
     --export=ALL,NSHARDS=$GRID_JOBS slurm/20_grid_run.sbatch)
j4=$(sbatch --parsable $EXTRA --cpus-per-task=$CPUS --dependency=afterok:$j3 --array=0-$((HARD_JOBS-1)) \
     --export=ALL,NSHARDS=$HARD_JOBS slurm/21_grid_hardchi.sbatch)
j5=$(sbatch --parsable $EXTRA --dependency=afterok:$j4 slurm/22_grid_analyze.sbatch)
j6=$(sbatch --parsable $EXTRA --dependency=afterany:$j5 slurm/30_timing.sbatch)
echo "submitted: confirm=$j1 confirm_analyze=$j2 grid=$j3 hardchi=$j4 grid_analyze=$j5 timing=$j6"
echo "查看: squeue -u $USER    日志: logs/"
