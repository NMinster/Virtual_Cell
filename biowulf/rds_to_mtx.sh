#!/bin/bash
# =============================================================================
# Virtual_Cell :: convert a Seurat / SCE .rds to an mtx folder, as a batch job
# =============================================================================
#     sbatch --output=logs/rds_%j.out biowulf/rds_to_mtx.sh <object.rds> <outdir> [assay] [layer]
#
# Uses Biowulf's R module (Seurat is installed there), not the conda env.
# `module` is a shell function on Biowulf and cannot be chained inside
# `sbatch --wrap="module load R && ..."` (Lmod eats the rest of the line), hence
# this script. Memory: a 6 GB Seurat object wants ~200 GB to be safe.
#SBATCH --cpus-per-task=4
#SBATCH --mem=240g
#SBATCH --time=6:00:00
#SBATCH --job-name=rds_to_mtx
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
module load R
R --version | head -1
Rscript R/rds_to_mtx.R "$@"
