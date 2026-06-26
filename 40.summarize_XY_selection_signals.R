#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
xdir <- file.path(base, "Genetics_Analysis/X_chromosome")

fst <- fread(file.path(xdir, "XY_FST_50kb_windows.tsv"))
pi  <- fread(file.path(xdir, "XY_pi_50kb_windows.tsv"))
td  <- fread(file.path(xdir, "XY_TajimaD_50kb_windows.tsv"))
het <- fread(file.path(xdir, "XY_Ho_50kb_windows.tsv"))

setnames(fst, c("CHR","window_start","window_end"), c("CHROM","BIN_START","BIN_END"))

x <- merge(
  fst[, .(CHROM, BIN_START, BIN_END, N_SNPs_FST = N_SNPs, Mean_FST, Max_FST, Top_SNP, Top_SNP_POS)],
  pi[, .(CHROM, BIN_START, BIN_END, PI_TT, PI_TX, delta_pi)],
  by = c("CHROM", "BIN_START", "BIN_END"),
  all = TRUE
)

x <- merge(
  x,
  td[, .(CHROM, BIN_START, BIN_END, TajimaD_TT, TajimaD_TX, delta_TajimaD)],
  by = c("CHROM", "BIN_START", "BIN_END"),
  all = TRUE
)

x <- merge(
  x,
  het[, .(CHROM, BIN_START, BIN_END, Ho_TT, Ho_TX, delta_Ho_TX_minus_TT)],
  by = c("CHROM", "BIN_START", "BIN_END"),
  all = TRUE
)

x[, Chromosome_type := fifelse(CHROM == "CM044047.1", "X",
                        fifelse(CHROM == "CM044048.1", "Y", "other"))]

x[, Evidence_score := 0L]
x[!is.na(Max_FST) & Max_FST >= 0.10, Evidence_score := Evidence_score + 1L]
x[!is.na(delta_pi) & delta_pi < 0, Evidence_score := Evidence_score + 1L]
x[!is.na(delta_TajimaD) & delta_TajimaD < 0, Evidence_score := Evidence_score + 1L]
x[!is.na(delta_Ho_TX_minus_TT) & delta_Ho_TX_minus_TT < 0, Evidence_score := Evidence_score + 1L]

setorder(x, -Evidence_score, -Max_FST)

out_all <- file.path(xdir, "XY_selection_signals_50kb_integrated.tsv")
out_top <- file.path(xdir, "XY_selection_signals_top50.tsv")

fwrite(x, out_all, sep = "\t", quote = FALSE, na = "NA")
fwrite(x[1:50], out_top, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("All:", out_all, "\n")
cat("Top:", out_top, "\n")
print(x[1:20, .(
  Chromosome_type, CHROM, BIN_START, BIN_END,
  Mean_FST, Max_FST, Top_SNP,
  delta_pi, delta_TajimaD, delta_Ho_TX_minus_TT,
  Evidence_score
)])
