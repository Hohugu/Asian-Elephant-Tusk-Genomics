#!/bin/bash
#SBATCH --job-name=gwas_plink
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gwas_plink.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gwas_plink.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

PRE="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS"
GWAS="/scratch/project_2000886/Hoedric/GWAS_2025/GWAS"

GENO="${PRE}/result/GWAS_qc_named"

PHENO="${PRE}/result/GWAS_pheno.txt"

COVAR5="${PRE}/result/GWAS_covar_PC5.txt"
COVAR10="${PRE}/result/GWAS_covar_PC10.txt"

OUTDIR="${GWAS}/result"

echo "[$(date)] GWAS logistic with PC1-PC5"

plink \
  --bfile "$GENO" \
  --allow-extra-chr \
  --allow-no-sex \
  --pheno "$PHENO" \
  --pheno-name tusk_num \
  --covar "$COVAR5" \
  --covar-name sex_bin,pop_bin,PC1,PC2,PC3,PC4,PC5 \
  --logistic hide-covar \
  --ci 0.95 \
  --out "${OUTDIR}/GWAS_logistic_PC5"

echo "[$(date)] GWAS logistic with PC1-PC10"

plink \
  --bfile "$GENO" \
  --allow-extra-chr \
  --allow-no-sex \
  --pheno "$PHENO" \
  --pheno-name tusk_num \
  --covar "$COVAR10" \
  --covar-name sex_bin,pop_bin,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --logistic hide-covar \
  --ci 0.95 \
  --out "${OUTDIR}/GWAS_logistic_PC10"

echo "[$(date)] PLINK logistic GWAS finished"
