#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fst_file <- file.path(base, "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_windows.tsv")
pi_file  <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_pi_50kb_comparison.tsv")
td_file  <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_comparison.tsv")
het_file <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv")
ld_file  <- file.path(base, "Genetics_Analysis/LD/tables/microLD_summary_TT_TX.tsv")

figdir <- file.path(base, "Genetics_Analysis/Diversity/figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

fst <- fread(fst_file)
pi  <- fread(pi_file)
td  <- fread(td_file)
het <- fread(het_file)
ld  <- fread(ld_file)

regions <- data.table(
  Region_ID = c("FST_region_18", "FST_region_9", "FST_region_56"),
  CHROM = c("CM044020.1", "CM044022.1", "CM044021.1"),
  START = c(77550001, 120450001, 136800001),
  END   = c(78100000, 121000000, 137350000)
)

for (i in seq_len(nrow(regions))) {

  r <- regions[i]
  chr <- r$CHROM
  start <- r$START
  end <- r$END
  rid <- r$Region_ID

  fst_sub <- fst[CHR == chr & window_start >= start & window_end <= end]
  pi_sub  <- pi[CHROM == chr & BIN_START >= start & BIN_END <= end]
  td_sub  <- td[CHROM == chr & BIN_START >= start & BIN_END <= end]
  het_sub <- het[CHROM == chr & BIN_START >= start & BIN_END <= end]

  png(
    file.path(figdir, paste0(rid, "_Campbell_like_zoom.png")),
    width = 1400,
    height = 1600
  )

  par(mfrow = c(5, 1), mar = c(3, 5, 3, 2))

  plot(
    fst_sub$window_start / 1e6,
    fst_sub$Mean_FST,
    pch = 20,
    type = "b",
    xlab = "",
    ylab = "Mean FST",
    main = paste0(rid, " | ", chr, ":", start, "-", end)
  )

  plot(
    pi_sub$BIN_START / 1e6,
    pi_sub$delta_pi,
    pch = 20,
    type = "b",
    xlab = "",
    ylab = expression(Delta*pi~"(TX-TT)")
  )
  abline(h = 0, lty = 2)

  plot(
    td_sub$BIN_START / 1e6,
    td_sub$delta_TajimaD,
    pch = 20,
    type = "b",
    xlab = "",
    ylab = expression(Delta~"TajimaD")
  )
  abline(h = 0, lty = 2)

  plot(
    het_sub$BIN_START / 1e6,
    het_sub$delta_Ho_TX_minus_TT,
    pch = 20,
    type = "b",
    xlab = "",
    ylab = expression(Delta~"Ho")
  )
  abline(h = 0, lty = 2)

  ld_row <- ld[region == rid]

  barplot(
    c(ld_row$mean_r2_TT, ld_row$mean_r2_TX),
    names.arg = c("TT", "TX"),
    ylab = expression(mean~r^2),
    main = "Local LD"
  )

  dev.off()

  cat("Wrote:", rid, "\n")
}

cat("Done\n")
