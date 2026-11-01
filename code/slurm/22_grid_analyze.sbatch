#!/bin/bash
#SBATCH --job-name=qecsb_grid_analyze
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=06:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_grid.py analyze --out results/paper1_grid --main-from results/paper1_v2/main_method.json
