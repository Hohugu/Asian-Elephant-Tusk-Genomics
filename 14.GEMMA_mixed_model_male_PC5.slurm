#!/bin/bash
#SBATCH --job-name=gemma_male_gwas
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_male_gwas.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/gemma_male_gwas.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"
PRE="${BASE}/PreGWAS/result"
GWAS="${BASE}/GWAS"

GENO="${PRE}/GWAS_qc_named"
PRUNED="${PRE}/GWAS_qc_named_pruned"

PHENO="${PRE}/GWAS_pheno.txt"
KEEP="${GWAS}/result/GEMMA_male_keep.txt"
COVAR="${GWAS}/result/GEMMA_male_covar_PC5.txt"

OUTDIR="${GWAS}/result/GEMMA_male"
mkdir -p "$OUTDIR"

echo "[$(date)] Create male-only full genotype"

plink \
  --bfile "$GENO" \
  --allow-extra-chr \
  --allow-no-sex \
  --keep "$KEEP" \
  --pheno "$PHENO" \
  --pheno-name tusk_num \
  --make-bed \
  --out "${OUTDIR}/GWAS_qc_named_male_pheno"

echo "[$(date)] Create male-only pruned genotype"

plink \
  --bfile "$PRUNED" \
  --allow-extra-chr \
  --allow-no-sex \
  --keep "$KEEP" \
  --pheno "$PHENO" \
  --pheno-name tusk_num \
  --make-bed \
  --out "${OUTDIR}/GWAS_qc_named_pruned_male_pheno"

cd "$OUTDIR"

echo "[$(date)] Compute male-only kinship matrix"

micromamba run -p /scratch/project_2000886/Hoedric/conda_envs/gemma_env \
  gemma \
  -bfile GWAS_qc_named_pruned_male_pheno \
  -gk 1 \
  -o tusk_male_kinship_pruned

echo "[$(date)] Run male-only GEMMA LMM"

micromamba run -p /scratch/project_2000886/Hoedric/conda_envs/gemma_env \
  gemma \
  -bfile GWAS_qc_named_male_pheno \
  -k output/tusk_male_kinship_pruned.cXX.txt \
  -c "$COVAR" \
  -lmm 4 \
  -o tusk_GEMMA_LMM_male_PC5

echo "[$(date)] Male-only GEMMA finished"
