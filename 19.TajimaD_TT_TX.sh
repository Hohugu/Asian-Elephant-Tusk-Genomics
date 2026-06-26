#!/bin/bash

set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
DIV=${BASE}/Genetics_Analysis/Diversity

VCF=${DIV}/GWAS_qc_named.vcf.gz

vcftools \
  --gzvcf ${VCF} \
  --keep ${DIV}/TT.vcf.samples \
  --TajimaD 50000 \
  --out ${DIV}/TT_TajimaD_50kb

vcftools \
  --gzvcf ${VCF} \
  --keep ${DIV}/TX.vcf.samples \
  --TajimaD 50000 \
  --out ${DIV}/TX_TajimaD_50kb
