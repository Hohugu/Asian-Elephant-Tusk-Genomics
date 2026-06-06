#!/bin/bash
#SBATCH --job-name=ld_pca
#SBATCH --output=ld_pca.out
#SBATCH --error=ld_pca.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail

module load plink/1.90b6.24

IN="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_qc_named"
OUTDIR="/scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result"

PRUNE="${OUTDIR}/GWAS_qc_named_pruned"
PCA="${OUTDIR}/GWAS_qc_named_pca"

echo "[$(date)] LD pruning"

plink \
  --bfile "$IN" \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.2 \
  --out "$PRUNE" \
  --memory 64000

echo "[$(date)] Create pruned dataset"

plink \
  --bfile "$IN" \
  --allow-extra-chr \
  --extract "${PRUNE}.prune.in" \
  --make-bed \
  --out "$PRUNE" \
  --memory 64000

echo "[$(date)] PCA"

plink \
  --bfile "$PRUNE" \
  --allow-extra-chr \
  --pca 20 \
  --out "$PCA" \
  --memory 64000

echo "[$(date)] LD pruning + PCA finished"
