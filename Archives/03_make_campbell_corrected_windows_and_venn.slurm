#!/bin/bash
#SBATCH --job-name=camp_corr_venn
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.err

set -euo pipefail

export BASE=/scratch/project_2000886/Hoedric/GWAS_2025
export CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction

export TMPDIR=$CAMP/Rtmp
export R_LIBS_USER=$CAMP/Rlibs

mkdir -p "$TMPDIR" "$R_LIBS_USER" "$CAMP/logs" "$CAMP/venn_campbell_corrected"

Rscript $CAMP/scripts/03_make_campbell_corrected_windows_and_venn.R
