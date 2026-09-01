#!/bin/bash
#SBATCH --job-name=cc_folds
#SBATCH --account=yitan
#SBATCH --partition=standard
#SBATCH --qos=part_qos_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --array=0-11%4
#SBATCH --output=/home/u10/hongbiaohuang/Log/cardiochrom_fold_%A_%a.out
#SBATCH --error=/home/u10/hongbiaohuang/Log/cardiochrom_fold_%A_%a.err

set -euo pipefail
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

PY=/xdisk/yitan/hongbiaohuang/software/venvs/virtual_cell_torch/bin/python
REPO=/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate

"$PY" -u "$REPO/scripts/export_frozen_bundles_on_puma.py" --fold-index "$SLURM_ARRAY_TASK_ID"
