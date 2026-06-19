#!/bin/bash
#SBATCH --job-name=fst_post
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/logs/fst_post.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/logs/fst_post.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2

set -euo pipefail

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"
TABLES="${BASE}/Genetics_Analysis/Selection/FST/tables"
LOGS="${BASE}/Genetics_Analysis/Selection/FST/logs"

export TMPDIR="${BASE}/tmp"
mkdir -p "$TMPDIR"

FST="${TABLES}/TT_vs_TX_genomewide.fst"

echo "[$(date)] Post-processing FST"

{
  head -1 "$FST"
  tail -n +2 "$FST" | sort -T "$TMPDIR" -S 20G -g -k5,5gr | head -100
} > "${TABLES}/TT_vs_TX_top100_FST_all.tsv"

{
  head -1 "$FST"
  tail -n +2 "$FST" | sort -T "$TMPDIR" -S 20G -g -k5,5gr | head -500
} > "${TABLES}/TT_vs_TX_top500_FST_all.tsv"

awk 'NR==1 || $1 ~ /^CM0440[2-4][0-9]\.1$/' "$FST" \
  > "${TABLES}/TT_vs_TX_genomewide_mainchr.fst"

{
  head -1 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst"
  tail -n +2 "${TABLES}/TT_vs_TX_genomewide_mainchr.fst" | sort -T "$TMPDIR" -S 20G -g -k5,5gr | head -100
} > "${TABLES}/TT_vs_TX_top100_FST_mainchr.tsv"

cat > "${LOGS}/make_FST_windows.R" << 'EOF'
library(data.table)

infile <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables/TT_vs_TX_genomewide.fst"
outdir <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables"

fst <- fread(infile)
fst_col <- names(fst)[grepl("FST", names(fst), ignore.case=TRUE)][1]

fst <- fst[!is.na(get(fst_col))]
fst <- fst[get(fst_col) >= 0]
fst[, POS := as.integer(tstrsplit(SNP, ":", fixed=TRUE)[[2]])]
fst <- fst[!is.na(POS)]

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

echo "[$(date)] Done"
