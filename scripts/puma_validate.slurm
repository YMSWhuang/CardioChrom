#!/bin/bash
#SBATCH --job-name=cc_validate
#SBATCH --account=yitan
#SBATCH --partition=standard
#SBATCH --qos=part_qos_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/u10/hongbiaohuang/Log/cardiochrom_validate_%j.out
#SBATCH --error=/home/u10/hongbiaohuang/Log/cardiochrom_validate_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

PY=/xdisk/yitan/hongbiaohuang/software/venvs/virtual_cell_torch/bin/python
REPO=/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate

"$PY" -u "$REPO/scripts/validate_against_scp3342_on_puma.py" --repo-root "$REPO" --n-jobs "$SLURM_CPUS_PER_TASK"
