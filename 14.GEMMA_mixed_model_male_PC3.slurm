#!/bin/bash
#SBATCH --job-name=gemma_male_PC3
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_male_PC3.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_male_PC3.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

export XDG_CACHE_HOME=/scratch/project_2000886/Hoedric/mamba_cache
export MAMBA_ROOT_PREFIX=/scratch/project_2000886/Hoedric/micromamba_root
export CONDA_PKGS_DIRS=/scratch/project_2000886/Hoedric/micromamba_root/pkgs

mkdir -p "$XDG_CACHE_HOME"
mkdir -p "$MAMBA_ROOT_PREFIX"
mkdir -p "$CONDA_PKGS_DIRS"

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"
GWAS="${BASE}/GWAS"

OUTDIR="${GWAS}/result/GEMMA_male"

COVAR="${GWAS}/result/GEMMA_male_covar_PC3.txt"

cd "$OUTDIR"

echo "[$(date)] Run male-only GEMMA PC3"

micromamba run \
-p /scratch/project_2000886/Hoedric/conda_envs/gemma_env \
gemma \
-bfile GWAS_qc_named_male_pheno \
-k output/tusk_male_kinship_pruned.cXX.txt \
-c "$COVAR" \
-lmm 4 \
-o tusk_GEMMA_LMM_male_PC3

echo "[$(date)] GEMMA male PC3 finished"

