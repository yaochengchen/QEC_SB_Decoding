#!/usr/bin/env bash
# 在已激活的 conda 环境里运行:  conda activate <你的环境>;  bash run.sh 1024
# 中断后再执行一次会接着跑（已完成的块会跳过）。
set -euo pipefail
cd "$(dirname "$0")"
WORKERS="${1:-1024}"
PY="${PYTHON:-python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "Python: $($PY -c 'import sys; print(sys.executable)')"

# 0) 依赖。ldpc 会顺带装 sinter，个别环境下 sinter 编译失败；失败时只装 ldpc 本体（脚本不需要 sinter）
$PY -m pip install -q -r requirements.txt
$PY -m pip install -q "ldpc>=2.3" || $PY -m pip install -q --no-deps "ldpc>=2.3" tqdm
mkdir -p results

# 1) 冒烟测试（几分钟）
$PY scripts/run_paper1_v2.py all --smoke 10 --d3-shots 200 --out results/smoke --workers "$WORKERS"
$PY scripts/run_paper1_grid.py all --scale 0.0005 --main CP_local_SOE --out results/smoke_grid --workers "$WORKERS"

# 2) 第一阶段：原确认集 12,000 shot（复现 + 基线 + 半径 2 消融），并按事先写好的规则选定主方法
$PY scripts/run_paper1_v2.py all --out results/paper1_v2 --workers "$WORKERS" 2>&1 | tee -a results/paper1_v2.log

# 闸门：复现不一致就停下，不跑网格
$PY - <<'PY'
import pandas as pd, sys
r = pd.read_csv("results/paper1_v2/reproduction_check.csv")
bad = r[(~r["match"].astype(bool)) | (r["r2_pools_rebuilt_identically"].astype(str) != "True")]
if len(bad):
    print("复现检查没有通过，先不要跑网格：\n", bad.to_string(index=False)); sys.exit(1)
print("复现检查通过。主方法:", open("results/paper1_v2/main_method.json").read())
PY

# 3) 第二阶段：p x d 网格（新 seed，主方法已固定）
$PY scripts/run_paper1_grid.py all --out results/paper1_grid --main-from results/paper1_v2/main_method.json \
    --workers "$WORKERS" 2>&1 | tee -a results/paper1_grid.log

echo "完成。计时请在没有其他任务的机器上单独运行："
echo "  $PY scripts/run_paper1_grid.py timing --out results/paper1_timing"
