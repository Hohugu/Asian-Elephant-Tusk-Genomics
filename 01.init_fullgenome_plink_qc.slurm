#!/bin/bash
#SBATCH --job-name=gwas_init_qc
#SBATCH --output=gwas_init_qc.out
#SBATCH --error=gwas_init_qc.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

VCF="/scratch/project_2000886/Hoedric/GWAS_2025/results/genotyped_vcf/all_contigs.raw.vcf.gz"
OUT="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_raw"

echo "[$(date)] Conversion VCF -> PLINK"

plink \
  --vcf "$VCF" \
  --make-bed \
  --allow-extra-chr \
  --vcf-half-call m \
  --memory 64000 \
  --out "$OUT"

echo "[$(date)] Allelic frequencies"

plink \
  --bfile "$OUT" \
  --freq \
  --allow-extra-chr \
  --out "${OUT}.freq_basic"

echo "[$(date)] Alleles counting"

plink \
  --bfile "$OUT" \
  --freq counts \
  --allow-extra-chr \
  --out "${OUT}.freq_counts"

echo "[$(date)] Genomics frequencies"

plink \
  --bfile "$OUT" \
  --freqx \
  --allow-extra-chr \
  --out "${OUT}.freqx"

echo "[$(date)] Missingness individus et variants"

plink \
  --bfile "$OUT" \
  --missing \
  --allow-extra-chr \
  --out "${OUT}.missing"

echo "[$(date)] Hardy-Weinberg"

plink \
  --bfile "$OUT" \
  --hardy \
  --allow-extra-chr \
  --out "${OUT}.hardy"

echo "[$(date)] Finished"
