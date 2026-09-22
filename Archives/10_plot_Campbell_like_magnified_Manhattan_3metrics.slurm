#!/bin/bash
#SBATCH --job-name=campbell_regplot
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/campbell_regplot.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/campbell_regplot.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r-env 2>/dev/null || true

CAMP=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction
mkdir -p "$CAMP/logs"

Rscript "$CAMP/scripts/10_plot_Campbell_like_magnified_Manhattan_3metrics.R"
