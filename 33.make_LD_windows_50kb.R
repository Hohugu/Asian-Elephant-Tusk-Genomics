#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

het_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/LD/regions/genomewide_50kb_windows_for_LD.tsv"
)

het <- fread(het_file)

win <- unique(
  het[
    N_SNPs >= 20,
    .(
      CHROM,
      BIN_START = as.integer(BIN_START),
      BIN_END = as.integer(BIN_END),
      N_SNPs
    )
  ]
)

setorder(win, CHROM, BIN_START)

fwrite(win, out_file, sep = "\t", quote = FALSE)

cat("Done\n")
cat("Windows:", nrow(win), "\n")
cat("Output:", out_file, "\n")
