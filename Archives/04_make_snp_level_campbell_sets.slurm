#!/bin/bash
#SBATCH --job-name=camp_snp_sets
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.err

set -euo pipefail

export BASE=/scratch/project_2000886/Hoedric/GWAS_2025
export CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction

mkdir -p $CAMP/logs
mkdir -p $CAMP/snp_level_campbell

python3 $CAMP/scripts/04_make_snp_level_campbell_sets.py \
  --camp $CAMP
