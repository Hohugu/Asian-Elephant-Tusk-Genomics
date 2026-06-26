#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fstpi_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_reduced_pi_overlapping_FST_annotated.tsv"
)

tajima_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_top100_more_negative_TX.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/FST_pi_TajimaD_convergence.tsv"
)

fstpi <- fread(fstpi_file)
taj <- fread(tajima_file)

fstpi[, key := paste(CHROM, START, END, sep=":")]
taj[, key := paste(CHROM, BIN_START + 1, BIN_END, sep=":")]

res <- merge(
  fstpi,
  taj[, .(
    key,
    TajimaD_TT,
    TajimaD_TX,
    delta_TajimaD
  )],
  by = "key",
  all = FALSE
)

setorder(res, delta_pi)

fwrite(res, out_file, sep="\t", quote=FALSE, na="NA")

cat("Done\n")
cat("Convergent regions:", uniqueN(res$Region_ID), "\n")
cat("Output:", out_file, "\n")
