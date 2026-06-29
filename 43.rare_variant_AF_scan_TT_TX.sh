#!/bin/bash
set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
DIV=$BASE/Genetics_Analysis/Diversity
OUT=$BASE/Genetics_Analysis/Functional_annotation/rare_variants

mkdir -p "$OUT"

VCF=$DIV/GWAS_qc_named.vcf.gz

vcftools --gzvcf "$VCF" \
  --keep "$DIV/TT.vcf.samples" \
  --freq2 \
  --out "$OUT/TT"

vcftools --gzvcf "$VCF" \
  --keep "$DIV/TX.vcf.samples" \
  --freq2 \
  --out "$OUT/TX"
