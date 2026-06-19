#!/bin/bash
#SBATCH --job-name=fst_TT_TX
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/logs/fst_TT_TX.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/logs/fst_TT_TX.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

set -euo pipefail

module load plink/1.90b6.24

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"
PRE="${BASE}/PreGWAS/result"
GA="${BASE}/Genetics_Analysis"

GENO="${PRE}/GWAS_qc_named"
CLUSTERS="${GA}/Population_structure/TT_TX_clusters.txt"

OUTDIR="${GA}/Selection/FST"
TABLES="${OUTDIR}/tables"
LOGS="${OUTDIR}/logs"

mkdir -p "$TABLES" "$LOGS"

echo "[$(date)] Starting genome-wide FST TT vs TX"

# 1. Individual SNP-level FST
plink \
  --bfile "$GENO" \
  --allow-extra-chr \
  --allow-no-sex \
  --within "$CLUSTERS" \
  --fst \
  --out "${TABLES}/TT_vs_TX_genomewide"

echo "[$(date)] SNP-level FST finished"

# 2. Top SNPs
{
  head -1 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst"
  tail -n +2 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst" | sort -g -k5,5nr | head -100
} > "${TABLES}/TT_vs_TX_top100_FST_mainchr.tsv"

{
  head -1 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst"
  tail -n +2 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst" | sort -g -k5,5nr | head -500
} > "${TABLES}/TT_vs_TX_top500_FST_mainchr.tsv"

# 3. Window summaries in R
cat > "${LOGS}/make_FST_windows.R" << 'EOF'
library(data.table)

infile <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables/TT_vs_TX_genomewide_mainchr.fst"
outdir <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables"

fst <- fread(infile)

# PLINK .fst columns usually: CHR SNP A1 A2 FST
# If your FST column has another name, detect it safely.
fst_col <- names(fst)[grepl("FST", names(fst), ignore.case=TRUE)][1]

fst <- fst[!is.na(get(fst_col))]
fst <- fst[get(fst_col) >= 0]

# Extract physical position from SNP ID: CM044020.1:85883031:A:T
fst[, POS := as.integer(tstrsplit(SNP, ":", fixed=TRUE)[[2]])]

make_windows <- function(dt, window_size, outfile) {
  x <- copy(dt)
  x[, window_start := floor((POS - 1) / window_size) * window_size + 1]
  x[, window_end := window_start + window_size - 1]

  res <- x[
    ,
    .(
      N_SNPs = .N,
      Mean_FST = mean(get(fst_col), na.rm=TRUE),
      Median_FST = median(get(fst_col), na.rm=TRUE),
      Max_FST = max(get(fst_col), na.rm=TRUE),
      Top_SNP = SNP[which.max(get(fst_col))]
    ),
    by=.(CHR, window_start, window_end)
  ]

  setorder(res, -Mean_FST)
  fwrite(res, file.path(outdir, outfile), sep="\t")
}

make_windows(fst, 10000,  "TT_vs_TX_FST_10kb_windows.tsv")
make_windows(fst, 50000,  "TT_vs_TX_FST_50kb_windows.tsv")
make_windows(fst, 100000, "TT_vs_TX_FST_100kb_windows.tsv")

cat("Done FST windows\n")
cat("SNPs used:", nrow(fst), "\n")
EOF

Rscript "${LOGS}/make_FST_windows.R"

echo "[$(date)] FST genome-wide scan finished"
