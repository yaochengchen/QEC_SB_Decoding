#!/bin/bash
# 过夜补充实验（在当前补充实验 10094 之后自动开始）：
#   AFTER=10094 bash slurm/submit_overnight.sh
#
# 1) 低 p 延伸：p=.04 全部 d（d<=9 40 万 shot，d>=11 100 万 shot），d=11/13 的 p=.05/.06（40 万 shot）
#    -> 所有方法逐 shot 配对，看次阈值区的 LER 随 d 的压低速度和有效距离
# 2) 新格点上的相关 BP+OSD；基准 CP-SB / CP-SA 的两组独立 solver seed 重复（d=9/11/13，p=.10/.14）
# 3) d>=11 难 shot 上 MPS chi=64 复核参照
# 4) 网格分析、补充分析全部重做（包含新格点）
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
source slurm/common.sh
DEP=""; [ -n "${AFTER:-}" ] && DEP="--dependency=afterany:$AFTER"
export QECSB_LOWP_EXT=1
G="python scripts/run_paper1_grid.py"; X="python scripts/run_paper1_extra.py"
A="sbatch --parsable $EXTRA --cpus-per-task=$CPUS --array=0-$((N-1)) --export=ALL,NSHARDS=$N"

j1=$(CMD="$G run --out results/paper1_grid" $A -J ext_grid $DEP slurm/generic_array.sbatch)
j2=$(CMD="$G hardchi --out results/paper1_grid" $A -J ext_hardchi --dependency=afterok:$j1 slurm/generic_array.sbatch)
j3=$(CMD="$X run --grid results/paper1_grid --out results/paper1_extra" $A -J ext_extra $DEP slurm/generic_array.sbatch)
j4=$(CMD="$X hard64 --grid results/paper1_grid --out results/paper1_extra" $A -J ext_hard64 --dependency=afterok:$j2 slurm/generic_array.sbatch)
j5=$(CMD="$G analyze --out results/paper1_grid --main-from results/paper1_v2/main_method.json" \
     sbatch --parsable $EXTRA --export=ALL -J ext_grid_analyze --dependency=afterok:$j2 slurm/generic_single.sbatch)
j6=$(CMD="$X analyze --grid results/paper1_grid --out results/paper1_extra" \
     sbatch --parsable $EXTRA --export=ALL -J ext_extra_analyze --dependency=afterok:$j3:$j4:$j5 slurm/generic_single.sbatch)
echo "submitted: grid=$j1 hardchi=$j2 extra=$j3 hard64=$j4 grid_analyze=$j5 extra_analyze=$j6"
