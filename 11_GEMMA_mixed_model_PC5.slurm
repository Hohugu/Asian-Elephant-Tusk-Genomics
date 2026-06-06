#!/bin/bash
#SBATCH --job-name=gemma_gwas
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_gwas.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_gwas.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

export MAMBA_ROOT_PREFIX=/scratch/project_2000886/Hoedric/micromamba_root
export CONDA_PKGS_DIRS=/scratch/project_2000886/Hoedric/micromamba_root/pkgs

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"

PRE="${BASE}/PreGWAS/result"
GWAS="${BASE}/GWAS"

GENO="${PRE}/GWAS_qc_named"
PRUNED="${PRE}/GWAS_qc_named_pruned"

PHENO="${PRE}/GWAS_pheno.txt"
COVAR="${GWAS}/result/GEMMA_covar_PC5.txt"

OUTDIR="${GWAS}/result/GEMMA"

mkdir -p "$OUTDIR"

echo "[$(date)] Create full genotype with phenotype"

plink --bfile "$GENO" --allow-extra-chr --allow-no-sex --pheno "$PHENO" --pheno-name tusk_num --make-bed --out "${OUTDIR}/GWAS_qc_named_pheno"

echo "[$(date)] Create pruned genotype with phenotype"

plink --bfile "$PRUNED" --allow-extra-chr --allow-no-sex --pheno "$PHENO" --pheno-name tusk_num --make-bed --out "${OUTDIR}/GWAS_qc_named_pruned_pheno"

cd "$OUTDIR"

echo "[$(date)] Compute kinship matrix"

micromamba run -p /scratch/project_2000886/Hoedric/conda_envs/gemma_env gemma -bfile GWAS_qc_named_pruned_pheno -gk 1 -o tusk_kinship_pruned

echo "[$(date)] Run GEMMA mixed model GWAS"

micromamba run -p /scratch/project_2000886/Hoedric/conda_envs/gemma_env gemma -bfile GWAS_qc_named_pheno -k output/tusk_kinship_pruned.cXX.txt -c "$COVAR" -lmm 4 -o tusk_GEMMA_LMM_PC5

echo "[$(date)] GEMMA GWAS finished"
