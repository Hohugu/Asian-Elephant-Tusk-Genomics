#!/bin/bash
#SBATCH --job-name=gwas_filter_male_qc
#SBATCH --output=gwas_filter_male_qc.out
#SBATCH --error=gwas_filter_male_qc.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

IN="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_raw_male"
OUT="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_male_qc_maf001_geno005_mind005"

plink \
  --bfile "$IN" \
  --allow-extra-chr \
  --snps-only just-acgt \
  --biallelic-only strict \
  --maf 0.01 \
  --geno 0.05 \
  --mind 0.05 \
  --make-bed \
  --memory 64000 \
  --out "$OUT"

plink --bfile "$OUT" --allow-extra-chr --missing --out "${OUT}.missing"
plink --bfile "$OUT" --allow-extra-chr --freq --out "${OUT}.freq"
plink --bfile "$OUT" --allow-extra-chr --het --out "${OUT}.het"

echo "[$(date)] QC filtering finished"
