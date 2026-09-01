#!/bin/bash
#SBATCH --job-name=cc_common
#SBATCH --account=yitan
#SBATCH --partition=standard
#SBATCH --qos=part_qos_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/u10/hongbiaohuang/Log/cardiochrom_common_%j.out
#SBATCH --error=/home/u10/hongbiaohuang/Log/cardiochrom_common_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

PY=/xdisk/yitan/hongbiaohuang/software/venvs/virtual_cell_torch/bin/python
REPO=/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate

"$PY" -u "$REPO/scripts/export_frozen_bundles_on_puma.py" --prepare-common
