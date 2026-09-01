#!/bin/bash
#SBATCH --job-name=cc_finalize
#SBATCH --account=yitan
#SBATCH --partition=standard
#SBATCH --qos=part_qos_standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/home/u10/hongbiaohuang/Log/cardiochrom_finalize_%j.out
#SBATCH --error=/home/u10/hongbiaohuang/Log/cardiochrom_finalize_%j.err

set -euo pipefail
export PYTHONUNBUFFERED=1

PY=/xdisk/yitan/hongbiaohuang/software/venvs/virtual_cell_torch/bin/python
REPO=/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate

"$PY" -u "$REPO/scripts/export_frozen_bundles_on_puma.py" --finalize
