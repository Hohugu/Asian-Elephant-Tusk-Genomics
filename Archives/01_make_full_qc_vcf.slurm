#!/bin/bash
#SBATCH --job-name=camp_vcf
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.err

set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction
BFILE=$BASE/GWAS/result/GEMMA_PC3/GWAS_qc_named_pheno
ENV=$CAMP/envs/raisd_build

export PATH=$ENV/bin:$CAMP/bin:$PATH
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-}

mkdir -p $CAMP/input $CAMP/logs

plink \
  --bfile "$BFILE" \
  --allow-extra-chr \
  --recode vcf bgz \
  --keep-allele-order \
  --out $CAMP/input/GWAS_qc_named_pheno

tabix -f -p vcf $CAMP/input/GWAS_qc_named_pheno.vcf.gz

bcftools query -l $CAMP/input/GWAS_qc_named_pheno.vcf.gz \
  > $CAMP/input/GWAS_qc_named_pheno.samples.txt

ls -lh $CAMP/input/GWAS_qc_named_pheno.vcf.gz*
