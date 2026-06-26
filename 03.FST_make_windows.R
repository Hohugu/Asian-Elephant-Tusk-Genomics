#!/usr/bin/env Rscript

library(data.table)

infile <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables/TT_vs_TX_genomewide.fst"
outdir <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables"

cat("Loading FST...\n")

fst <- fread(infile)

cat("Columns:\n")
print(names(fst))

fst <- fst[!is.na(FST)]
fst <- fst[FST >= 0]

make_windows <- function(dt, window_size, outfile) {

  cat("Making", window_size, "bp windows...\n")

  x <- copy(dt)

  x[, window_start := floor((POS - 1) / window_size) * window_size + 1]
  x[, window_end := window_start + window_size - 1]

  res <- x[
    ,
    .(
      N_SNPs = .N,
      Mean_FST = mean(FST, na.rm=TRUE),
      Median_FST = median(FST, na.rm=TRUE),
      Max_FST = max(FST, na.rm=TRUE),
      Top_SNP = SNP[which.max(FST)],
      Top_SNP_POS = POS[which.max(FST)]
    ),
    by=.(CHR, window_start, window_end)
  ]

  setorder(res, -Mean_FST)

  fwrite(res, file.path(outdir, outfile), sep="\t")

  cat("Wrote:", file.path(outdir, outfile), "\n")
}

make_windows(fst, 10000,  "TT_vs_TX_FST_10kb_windows.tsv")
make_windows(fst, 50000,  "TT_vs_TX_FST_50kb_windows.tsv")
make_windows(fst, 100000, "TT_vs_TX_FST_100kb_windows.tsv")

cat("Done\n")
cat("SNPs used:", nrow(fst), "\n")
