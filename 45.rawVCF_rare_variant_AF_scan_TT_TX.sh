#!/bin/bash
set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
DIV=$BASE/Genetics_Analysis/Diversity
OUT=$BASE/Genetics_Analysis/Functional_annotation/rare_variants_rawVCF

mkdir -p "$OUT"

VCF=$BASE/results/genotyped_vcf/all_contigs.raw.vcf.gz

echo "Checking keep files..."
echo "TT overlap:"
comm -12 <(bcftools query -l "$VCF" | sort) <(sort "$DIV/TT.raw.keep") | wc -l
echo "TX overlap:"
comm -12 <(bcftools query -l "$VCF" | sort) <(sort "$DIV/TX.raw.keep") | wc -l

vcftools --gzvcf "$VCF" \
  --keep "$DIV/TT.raw.keep" \
  --freq2 \
  --out "$OUT/TT_raw"

vcftools --gzvcf "$VCF" \
  --keep "$DIV/TX.raw.keep" \
  --freq2 \
  --out "$OUT/TX_raw"

echo "Done"
echo "Output: $OUT"
