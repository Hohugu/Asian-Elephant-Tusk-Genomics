#!/bin/bash
#SBATCH --job-name=heterozygosity
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load plink

BASE=/scratch/project_2000886/Hoedric/GWAS_2025

plink \
  --bfile ${BASE}/PreGWAS/result/GWAS_qc_named \
  --allow-extra-chr \
  --allow-no-sex \
  --het \
  --out ${BASE}/Genetics_Analysis/Diversity/TT_TX_heterozygosity
