#!/bin/bash
# 在登录节点上执行一次（需要能联网）：  bash slurm/setup_env.sh
set -euo pipefail
cd "$(dirname "$0")/.."
module load Core/anaconda3-2024.02-1
eval "$(conda shell.bash hook)"
ENV=qecsb
if ! conda env list | grep -q "^$ENV "; then
    conda create -y -n "$ENV" python=3.11
    # home 空间不够时改用:  conda create -y -p /你的大盘/envs/qecsb python=3.11 ，并把 env.sh 里的 ENV 改成该路径
fi
conda activate "$ENV"
python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install "ldpc>=2.3" || python -m pip install --no-deps "ldpc>=2.3" tqdm   # sinter 装不上时的退路
python - <<'PY'
import numpy, scipy, pandas, networkx, stim, pymatching, qecsim, matplotlib, tabulate
import sys; sys.path.insert(0, "scripts")
from run_paper1_v2 import import_bposd; import_bposd()
print("环境 OK:", sys.executable)
PY
