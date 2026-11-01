#!/bin/bash
#SBATCH --job-name=qecsb_extra_analyze
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --exclude=n013
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
source slurm/env.sh
set -euo pipefail
python scripts/run_paper1_extra.py analyze --grid results/paper1_grid --out results/paper1_extra
