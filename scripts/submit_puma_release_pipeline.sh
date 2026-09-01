#!/bin/bash
set -euo pipefail

REPO=/xdisk/yitan/hongbiaohuang/CardioChrom_release_candidate
mkdir -p /home/u10/hongbiaohuang/Log

PREP_JOB=$(sbatch --parsable "$REPO/scripts/puma_export_common.slurm")
FOLD_JOB=$(sbatch --parsable --dependency="afterok:${PREP_JOB}" "$REPO/scripts/puma_export_folds.slurm")
FINAL_JOB=$(sbatch --parsable --dependency="afterok:${FOLD_JOB}" "$REPO/scripts/puma_export_finalize.slurm")
VALIDATE_JOB=$(sbatch --parsable --dependency="afterok:${FINAL_JOB}" "$REPO/scripts/puma_validate.slurm")

echo "prepare_common=${PREP_JOB}"
echo "export_folds=${FOLD_JOB}"
echo "finalize=${FINAL_JOB}"
echo "validate=${VALIDATE_JOB}"
