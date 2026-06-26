#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
ld_file <- file.path(base, "Genetics_Analysis/LD/tables/microLD_summary_TT_TX.tsv")
figdir <- file.path(base, "Genetics_Analysis/LD/figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

dat <- fread(ld_file)

dat[, region := factor(region, levels = region[order(delta_mean_r2_TX_minus_TT)])]

png(file.path(figdir, "microLD_mean_r2_TT_TX.png"), width = 1400, height = 900)

par(mar = c(10, 5, 4, 2))

bp <- barplot(
  rbind(dat$mean_r2_TT, dat$mean_r2_TX),
  beside = TRUE,
  names.arg = dat$region,
  las = 2,
  ylab = expression(mean~r^2),
  main = "Local LD comparison between TT and TX",
  ylim = c(0, max(dat$mean_r2_TT, dat$mean_r2_TX) * 1.2)
)

legend(
  "topright",
  legend = c("TT", "TX"),
  fill = gray.colors(2),
  bty = "n"
)

dev.off()

png(file.path(figdir, "microLD_delta_mean_r2_TX_minus_TT.png"), width = 1400, height = 900)

par(mar = c(10, 5, 4, 2))

barplot(
  dat$delta_mean_r2_TX_minus_TT,
  names.arg = dat$region,
  las = 2,
  ylab = expression(Delta~mean~r^2~"(TX - TT)"),
  main = "Difference in local LD between TX and TT"
)

abline(h = 0, lty = 2)

dev.off()

png(file.path(figdir, "microLD_highLD_pairs_TT_TX.png"), width = 1400, height = 900)

par(mar = c(10, 5, 4, 2))

barplot(
  rbind(dat$n_r2_over_0.8_TT, dat$n_r2_over_0.8_TX),
  beside = TRUE,
  names.arg = dat$region,
  las = 2,
  ylab = "Number of SNP pairs with r² >= 0.8",
  main = "High-LD SNP pairs by region"
)

legend(
  "topright",
  legend = c("TT", "TX"),
  fill = gray.colors(2),
  bty = "n"
)

dev.off()

cat("Done\n")
cat("Figures written to:", figdir, "\n")
