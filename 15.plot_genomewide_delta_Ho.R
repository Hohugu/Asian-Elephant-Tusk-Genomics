#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
div <- file.path(base, "Genetics_Analysis/Diversity")
figdir <- file.path(div, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

het <- fread(file.path(div, "TT_TX_windowed_observed_heterozygosity_50kb.tsv"))

het <- het[N_SNPs >= 10]
het[, BIN_START := as.numeric(BIN_START)]
het[, BIN_END := as.numeric(BIN_END)]

het[, chr_order := match(CHROM, unique(CHROM))]
setorder(het, chr_order, BIN_START)

chr_sizes <- het[, .(chr_len = max(BIN_END)), by = .(CHROM, chr_order)]
setorder(chr_sizes, chr_order)
chr_sizes[, offset := c(0, cumsum(as.numeric(chr_len))[-.N])]

het <- merge(het, chr_sizes[, .(CHROM, offset)], by = "CHROM", all.x = TRUE)
het[, BPcum := BIN_START + offset]

axis_df <- het[
  ,
  .(center = (min(BPcum) + max(BPcum)) / 2),
  by = .(CHROM, chr_order)
]
setorder(axis_df, chr_order)

axis_df[, label := fifelse(
  grepl("^CM0440", CHROM),
  sub("^CM0440", "", sub("\\.1$", "", CHROM)),
  ""
)]

p01_low <- quantile(het$delta_Ho_TX_minus_TT, 0.001, na.rm = TRUE)
p1_low  <- quantile(het$delta_Ho_TX_minus_TT, 0.01, na.rm = TRUE)

fwrite(
  data.table(
    threshold = c("bottom_0.1_percent", "bottom_1_percent"),
    delta_Ho = c(p01_low, p1_low)
  ),
  file.path(div, "TT_TX_delta_Ho_50kb_thresholds.tsv"),
  sep = "\t"
)

top100 <- het[order(delta_Ho_TX_minus_TT)][1:100]
fwrite(
  top100,
  file.path(div, "TT_TX_delta_Ho_50kb_top100_reduced_TX.tsv"),
  sep = "\t",
  quote = FALSE
)

png(
  file.path(figdir, "TT_TX_delta_Ho_50kb_genomewide.png"),
  width = 2400,
  height = 900
)

plot(
  het$BPcum,
  het$delta_Ho_TX_minus_TT,
  pch = 20,
  cex = 0.35,
  col = ifelse(het$delta_Ho_TX_minus_TT <= p01_low, "blue",
               ifelse(het$delta_Ho_TX_minus_TT <= p1_low, "orange", "grey70")),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(Delta~Ho~"(TX - TT)"),
  main = "Genome-wide observed heterozygosity difference, 50 kb windows"
)

axis(1, at = axis_df$center, labels = axis_df$label, las = 2, cex.axis = 0.8)
abline(h = 0, lty = 2)
abline(h = p1_low, lty = 2, col = "orange")
abline(h = p01_low, lty = 2, col = "blue")

legend(
  "bottomright",
  legend = c("bottom 1%", "bottom 0.1%"),
  col = c("orange", "blue"),
  pch = 20,
  bty = "n"
)

dev.off()

cat("Done\n")
cat("Windows:", nrow(het), "\n")
cat("Bottom 1%:", p1_low, "\n")
cat("Bottom 0.1%:", p01_low, "\n")
cat("Figure:", file.path(figdir, "TT_TX_delta_Ho_50kb_genomewide.png"), "\n")
