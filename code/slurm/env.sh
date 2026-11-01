# 被各个作业脚本 source；不要单独运行
module load Core/anaconda3-2024.02-1
eval "$(conda shell.bash hook)"
conda activate qecsb
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1
cd "$SLURM_SUBMIT_DIR" || { echo "Failed to cd to $SLURM_SUBMIT_DIR"; exit 1; }
W=${SLURM_CPUS_PER_TASK:-1}
