#!/bin/bash
#SBATCH --job-name=qecsb_smoke
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --time=01:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_v2.py all --smoke 10 --d3-shots 200 --out results/smoke --workers "$W"
python scripts/run_paper1_grid.py all --scale 0.0005 --main CP_local_SOE --out results/smoke_grid --workers "$W"
python scripts/run_paper1_grid.py estimate --out results/estimate --workers 2048 | tee results/estimate.txt
echo SMOKE_OK
