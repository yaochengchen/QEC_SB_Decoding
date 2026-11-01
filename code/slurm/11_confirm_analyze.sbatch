#!/bin/bash
# 阶段 1 分析 + 复现闸门。闸门不通过时本作业以失败退出，后面依赖它的网格作业不会启动。
#SBATCH --job-name=qecsb_confirm_analyze
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_v2.py check  --out "${OUT:-results/paper1_v2}" --base-seed "${BASE_SEED:-2027101300}"
python scripts/run_paper1_v2.py mps-d3 --out "${OUT:-results/paper1_v2}" --base-seed "${BASE_SEED:-2027101300}"
python scripts/run_paper1_v2.py analyze --out "${OUT:-results/paper1_v2}" --base-seed "${BASE_SEED:-2027101300}"
python - <<'PY'
import pandas as pd, sys
r = pd.read_csv("" + __import__("os").environ.get("OUT","results/paper1_v2") + "/reproduction_check.csv")
bad = r[r["r2_pools_rebuilt_identically"].astype(str) != "True"]   # 内部一致性检查；不与旧论文表比较
if len(bad):
    print("CONSISTENCY CHECK FAILED, grid will not start:\n", bad.to_string(index=False)); sys.exit(1)
print("consistency OK; main method:", open(__import__("os").environ.get("OUT","results/paper1_v2") + "/main_method.json").read())
PY
