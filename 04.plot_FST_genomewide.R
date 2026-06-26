#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST"

figdir <- file.path(base, "figures")
tabdir <- file.path(base, "tables")

dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

plot_fst <- function(infile, outfile, title) {

  cat("Loading:", infile, "\n")

  fst <- fread(file.path(tabdir, infile))

  fst <- fst[!is.na(Mean_FST)]
  fst <- fst[N_SNPs >= 2]

  p99  <- quantile(fst$Mean_FST, 0.99,  na.rm = TRUE)
  p999 <- quantile(fst$Mean_FST, 0.999, na.rm = TRUE)

  cat("99th percentile:", p99, "\n")
  cat("99.9th percentile:", p999, "\n")

  fst[, category := "background"]
  fst[Mean_FST >= p99,  category := "p99"]
  fst[Mean_FST >= p999, category := "p999"]

  cols <- c(
    background = "grey75",
    p99 = "orange",
    p999 = "blue"
  )

  sizes <- ifelse(
    fst$category == "p999", 0.9,
    ifelse(fst$category == "p99", 0.55, 0.22)
  )

  fst[, chr_type := fifelse(grepl("^CM", CHR), "CM", "scaffold")]

  chr_order <- unique(fst[order(chr_type, CHR)]$CHR)
  fst[, CHR_index := match(CHR, chr_order)]

  setorder(fst, CHR_index, window_start)

  chr_sizes <- fst[
    ,
    .(chr_len = max(window_end, na.rm = TRUE)),
    by = .(CHR, CHR_index)
  ]

  setorder(chr_sizes, CHR_index)

  chr_sizes[, offset := c(0, cumsum(chr_len)[-.N])]

  fst <- merge(
    fst,
    chr_sizes[, .(CHR, offset)],
    by = "CHR",
    all.x = TRUE
  )

  fst[, BPcum := window_start + offset]

  axis_df <- fst[
    ,
    .(center = (min(BPcum) + max(BPcum)) / 2),
    by = .(CHR, CHR_index)
  ]

  setorder(axis_df, CHR_index)

  axis_df[, label := fifelse(
    grepl("^CM0440", CHR),
    sub("^CM0440", "", sub("\\.1$", "", CHR)),
    ""
  )]

  chr_boundaries <- chr_sizes$offset[-1]

  png(file.path(figdir, outfile), width = 2400, height = 950)

  par(mar = c(6, 5, 4, 2))

  plot(
    fst$BPcum,
    fst$Mean_FST,
    pch = 20,
    cex = sizes,
    col = cols[fst$category],
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression("Mean " * F[ST]),
    main = title
  )

  abline(v = chr_boundaries, col = "grey90", lwd = 0.5)

  axis(
    1,
    at = axis_df$center,
    labels = axis_df$label,
    las = 2,
    cex.axis = 0.8
  )

  abline(h = p99,  col = "orange", lty = 2, lwd = 2)
  abline(h = p999, col = "blue",   lty = 2, lwd = 2)

  legend(
    "topright",
    legend = c("Background", "Top 1%", "Top 0.1%", "99th percentile", "99.9th percentile"),
    col = c("grey75", "orange", "blue", "orange", "blue"),
    pch = c(20, 20, 20, NA, NA),
    lty = c(NA, NA, NA, 2, 2),
    lwd = c(NA, NA, NA, 2, 2),
    pt.cex = c(0.7, 0.9, 1.1, NA, NA),
    bty = "n"
  )

  dev.off()

  thresholds <- data.table(
    Window_file = infile,
    Threshold = c("p99", "p999"),
    Mean_FST = c(p99, p999)
  )

  fwrite(
    thresholds,
    file.path(tabdir, sub(".png", "_thresholds.tsv", outfile)),
    sep = "\t"
  )

  candidate_windows <- fst[
    Mean_FST >= p999,
    .(
      CHR,
      window_start,
      window_end,
      N_SNPs,
      Mean_FST,
      Median_FST,
      Max_FST,
      Top_SNP,
      Top_SNP_POS
    )
  ]

  setorder(candidate_windows, -Mean_FST)

  fwrite(
    candidate_windows,
    file.path(tabdir, sub(".png", "_candidate_windows_p999.tsv", outfile)),
    sep = "\t"
  )

  cat("Wrote figure:", file.path(figdir, outfile), "\n")
  cat("Wrote candidate windows:", file.path(tabdir, sub(".png", "_candidate_windows_p999.tsv", outfile)), "\n")
}

plot_fst(
  infile = "TT_vs_TX_FST_10kb_windows.tsv",
  outfile = "TT_vs_TX_FST_10kb_genomewide_highlight.png",
  title = "Genome-wide FST scan TT vs TX, 10 kb windows"
)

plot_fst(
  infile = "TT_vs_TX_FST_50kb_windows.tsv",
  outfile = "TT_vs_TX_FST_50kb_genomewide_highlight.png",
  title = "Genome-wide FST scan TT vs TX, 50 kb windows"
)

cat("Done\n")
