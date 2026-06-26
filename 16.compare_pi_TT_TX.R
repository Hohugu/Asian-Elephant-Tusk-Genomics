#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Diversity"

tt_file <- file.path(base, "TT_pi_50kb.windowed.pi")
tx_file <- file.path(base, "TX_pi_50kb.windowed.pi")

figdir <- file.path(base, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

tt <- fread(tt_file)
tx <- fread(tx_file)

setnames(tt, "PI", "PI_TT")
setnames(tx, "PI", "PI_TX")

pi <- merge(
  tt[, .(CHROM, BIN_START, BIN_END, N_VARIANTS_TT = N_VARIANTS, PI_TT)],
  tx[, .(CHROM, BIN_START, BIN_END, N_VARIANTS_TX = N_VARIANTS, PI_TX)],
  by = c("CHROM", "BIN_START", "BIN_END")
)

pi <- pi[!is.na(PI_TT) & !is.na(PI_TX)]
pi <- pi[N_VARIANTS_TT >= 10 & N_VARIANTS_TX >= 10]

pi[, delta_pi := PI_TX - PI_TT]
pi[, pi_ratio_TX_TT := PI_TX / PI_TT]
pi[, log2_pi_ratio_TX_TT := log2(pi_ratio_TX_TT)]

pi[!is.finite(log2_pi_ratio_TX_TT), log2_pi_ratio_TX_TT := NA_real_]

summary <- data.table(
  metric = c(
    "mean_PI_TT",
    "mean_PI_TX",
    "median_PI_TT",
    "median_PI_TX",
    "mean_delta_PI_TX_minus_TT",
    "median_delta_PI_TX_minus_TT",
    "wilcox_p_PI_TT_vs_TX"
  ),
  value = c(
    mean(pi$PI_TT, na.rm = TRUE),
    mean(pi$PI_TX, na.rm = TRUE),
    median(pi$PI_TT, na.rm = TRUE),
    median(pi$PI_TX, na.rm = TRUE),
    mean(pi$delta_pi, na.rm = TRUE),
    median(pi$delta_pi, na.rm = TRUE),
    wilcox.test(pi$PI_TT, pi$PI_TX, paired = TRUE)$p.value
  )
)

fwrite(pi, file.path(base, "TT_TX_pi_50kb_comparison.tsv"), sep = "\t")
fwrite(summary, file.path(base, "TT_TX_pi_50kb_summary.tsv"), sep = "\t")

# Top windows with reduced diversity in TX
reduced_TX <- pi[order(delta_pi)][1:100]
fwrite(reduced_TX, file.path(base, "TT_TX_pi_50kb_top100_reduced_TX.tsv"), sep = "\t")

# Top windows with reduced diversity in TT
reduced_TT <- pi[order(-delta_pi)][1:100]
fwrite(reduced_TT, file.path(base, "TT_TX_pi_50kb_top100_reduced_TT.tsv"), sep = "\t")

# Genome-wide plot of delta pi
pi[, chr_order := match(CHROM, unique(CHROM))]
setorder(pi, chr_order, BIN_START)

chr_sizes <- pi[, .(chr_len = max(BIN_END)), by = .(CHROM, chr_order)]
setorder(chr_sizes, chr_order)
chr_sizes[, offset := c(0, cumsum(chr_len)[-.N])]

pi <- merge(pi, chr_sizes[, .(CHROM, offset)], by = "CHROM", all.x = TRUE)
pi[, BPcum := BIN_START + offset]

axis_df <- pi[
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

png(file.path(figdir, "TT_TX_delta_pi_50kb_genomewide.png"), width = 2400, height = 900)

plot(
  pi$BPcum,
  pi$delta_pi,
  pch = 20,
  cex = 0.3,
  col = ifelse(pi$delta_pi < 0, "blue", "orange"),
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(pi[TX] - pi[TT]),
  main = "Genome-wide nucleotide diversity difference, 50 kb windows"
)

axis(1, at = axis_df$center, labels = axis_df$label, las = 2, cex.axis = 0.8)
abline(h = 0, lty = 2, lwd = 2)
dev.off()

# Boxplot PI
png(file.path(figdir, "TT_TX_pi_50kb_boxplot.png"), width = 900, height = 900)

boxplot(
  list(TT = pi$PI_TT, TX = pi$PI_TX),
  ylab = expression(pi),
  main = "Nucleotide diversity, 50 kb windows",
  outline = FALSE
)

dev.off()

cat("Done\n")
cat("Windows retained:", nrow(pi), "\n")
cat("Output:", file.path(base, "TT_TX_pi_50kb_comparison.tsv"), "\n")
print(summary)
