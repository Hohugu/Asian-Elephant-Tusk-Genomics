#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Diversity"
figdir <- file.path(base, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

tt <- fread(file.path(base, "TT_TajimaD_50kb.Tajima.D"))
tx <- fread(file.path(base, "TX_TajimaD_50kb.Tajima.D"))

setnames(tt, "TajimaD", "TajimaD_TT")
setnames(tx, "TajimaD", "TajimaD_TX")

dat <- merge(
  tt[, .(CHROM, BIN_START, N_SNPS_TT = N_SNPS, TajimaD_TT)],
  tx[, .(CHROM, BIN_START, N_SNPS_TX = N_SNPS, TajimaD_TX)],
  by = c("CHROM", "BIN_START")
)

dat <- dat[N_SNPS_TT >= 10 & N_SNPS_TX >= 10]

dat[, delta_TajimaD := TajimaD_TX - TajimaD_TT]
dat[, BIN_END := BIN_START + 50000]

summary <- data.table(
  metric = c(
    "mean_TajimaD_TT",
    "mean_TajimaD_TX",
    "median_TajimaD_TT",
    "median_TajimaD_TX",
    "mean_delta_TX_minus_TT",
    "median_delta_TX_minus_TT",
    "wilcox_p_TT_vs_TX"
  ),
  value = c(
    mean(dat$TajimaD_TT, na.rm = TRUE),
    mean(dat$TajimaD_TX, na.rm = TRUE),
    median(dat$TajimaD_TT, na.rm = TRUE),
    median(dat$TajimaD_TX, na.rm = TRUE),
    mean(dat$delta_TajimaD, na.rm = TRUE),
    median(dat$delta_TajimaD, na.rm = TRUE),
    wilcox.test(dat$TajimaD_TT, dat$TajimaD_TX, paired = TRUE)$p.value
  )
)

fwrite(dat, file.path(base, "TT_TX_TajimaD_50kb_comparison.tsv"), sep = "\t")
fwrite(summary, file.path(base, "TT_TX_TajimaD_50kb_summary.tsv"), sep = "\t")

fwrite(dat[order(delta_TajimaD)][1:100],
       file.path(base, "TT_TX_TajimaD_50kb_top100_more_negative_TX.tsv"),
       sep = "\t")

fwrite(dat[order(-delta_TajimaD)][1:100],
       file.path(base, "TT_TX_TajimaD_50kb_top100_more_positive_TX.tsv"),
       sep = "\t")

dat[, chr_order := match(CHROM, unique(CHROM))]
setorder(dat, chr_order, BIN_START)

chr_sizes <- dat[, .(chr_len = max(BIN_END)), by = .(CHROM, chr_order)]
setorder(chr_sizes, chr_order)
chr_sizes[, offset := c(0, cumsum(as.numeric(chr_len))[-.N])]

dat <- merge(dat, chr_sizes[, .(CHROM, offset)], by = "CHROM", all.x = TRUE)
dat[, BPcum := BIN_START + offset]

axis_df <- dat[
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

png(file.path(figdir, "TT_TX_delta_TajimaD_50kb_genomewide.png"),
    width = 2400, height = 900)

plot(
  dat$BPcum,
  dat$delta_TajimaD,
  pch = 20,
  cex = 0.3,
  col = ifelse(dat$delta_TajimaD < 0, "blue", "orange"),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(Delta~"Tajima's D (TX - TT)"),
  main = "Genome-wide Tajima's D difference, 50 kb windows"
)

axis(1, at = axis_df$center, labels = axis_df$label, las = 2, cex.axis = 0.8)
abline(h = 0, lty = 2, lwd = 2)

dev.off()

png(file.path(figdir, "TT_TX_TajimaD_50kb_boxplot.png"),
    width = 900, height = 900)

boxplot(
  list(TT = dat$TajimaD_TT, TX = dat$TajimaD_TX),
  ylab = "Tajima's D",
  main = "Tajima's D, 50 kb windows",
  outline = FALSE
)

dev.off()

cat("Done\n")
cat("Windows retained:", nrow(dat), "\n")
print(summary)
