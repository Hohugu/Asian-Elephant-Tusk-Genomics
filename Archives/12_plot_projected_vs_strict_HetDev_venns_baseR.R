BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
INDIR <- file.path(CAMP, "strict_HetDev_venn")
OUTDIR <- file.path(CAMP, "strict_HetDev_venn_figures_baseR")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

ellipse <- function(cx, cy, rx, ry, angle = 0, n = 300) {
  t <- seq(0, 2*pi, length.out = n)
  a <- angle * pi / 180
  x <- rx * cos(t)
  y <- ry * sin(t)
  data.frame(
    x = cx + x * cos(a) - y * sin(a),
    y = cy + x * sin(a) + y * cos(a)
  )
}

read_counts <- function(label, mode) {
  f <- file.path(INDIR, paste0("pattern_counts_", label, "_", mode, ".tsv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  x <- read.delim(f, stringsAsFactors = FALSE)
  x <- x[x$N > 0, ]
  x
}

pattern_positions <- data.frame(
  pattern = c(
    "FST", "DXY", "muLD", "HetDev_TX",
    "FST&DXY", "FST&muLD", "FST&HetDev_TX",
    "DXY&muLD", "DXY&HetDev_TX", "muLD&HetDev_TX",
    "FST&DXY&muLD", "FST&DXY&HetDev_TX",
    "FST&muLD&HetDev_TX", "DXY&muLD&HetDev_TX",
    "FST&DXY&muLD&HetDev_TX"
  ),
  x = c(
    -1.75,  1.75, -1.75,  1.75,
     0.00, -1.45,  0.00,
     0.00,  1.45,  0.00,
    -0.45,  0.45,
    -0.45,  0.45,
     0.00
  ),
  y = c(
     1.10,  1.10, -1.10, -1.10,
     1.25,  0.00,  0.35,
    -0.35,  0.00, -1.25,
     0.45,  0.45,
    -0.45, -0.45,
     0.00
  ),
  stringsAsFactors = FALSE
)

plot_venn_like <- function(counts, label, mode, ext = "png") {
  outfile <- file.path(
    OUTDIR,
    paste0("Campbell_SNP_level_Venn_like_", label, "_", mode, ".", ext)
  )

  if (ext == "png") {
    png(outfile, width = 2200, height = 1800, res = 250)
  } else {
    pdf(outfile, width = 8.8, height = 7.2)
  }

  par(mar = c(1, 1, 4, 1))
  plot(
    NA,
    xlim = c(-2.7, 2.7),
    ylim = c(-1.9, 1.9),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = paste0("Campbell SNP-level Venn-like: ", label, " ", mode)
  )

  sets <- list(
    list(name = "FST",       cx = -0.75, cy =  0.45, rx = 1.25, ry = 0.85, angle =  30, col = "lightblue"),
    list(name = "DXY",       cx =  0.75, cy =  0.45, rx = 1.25, ry = 0.85, angle = -30, col = "moccasin"),
    list(name = "muLD",      cx = -0.75, cy = -0.45, rx = 1.25, ry = 0.85, angle = -30, col = "palegreen"),
    list(name = "HetDev_TX", cx =  0.75, cy = -0.45, rx = 1.25, ry = 0.85, angle =  30, col = "mistyrose")
  )

  for (s in sets) {
    e <- ellipse(s$cx, s$cy, s$rx, s$ry, s$angle)
    polygon(e$x, e$y, col = adjustcolor(s$col, alpha.f = 0.45), border = "grey30", lwd = 1.4)
  }

  text(-2.15,  1.45, "FST", font = 2, cex = 1.1)
  text( 2.15,  1.45, "DXY", font = 2, cex = 1.1)
  text(-2.15, -1.45, "muLD", font = 2, cex = 1.1)
  text( 2.15, -1.45, "HetDev_TX", font = 2, cex = 1.1)

  lab <- merge(counts, pattern_positions, by = "pattern", all.x = TRUE)
  lab <- lab[!is.na(lab$x), ]

  for (i in seq_len(nrow(lab))) {
    txt <- paste0(lab$pattern[i], "\n", format(lab$N[i], big.mark = ",", scientific = FALSE))
    text(lab$x[i], lab$y[i], txt, cex = 0.72)
  }

  mtext("Non-proportional schematic; labels are exclusive SNP counts", side = 3, line = 0.4, cex = 0.85)

  dev.off()
  message("Wrote: ", outfile)
}

plot_upset_bar <- function(counts, label, mode, ext = "png") {
  outfile <- file.path(
    OUTDIR,
    paste0("Campbell_SNP_level_UpSet_bar_", label, "_", mode, ".", ext)
  )

  counts <- counts[order(counts$n_metrics, counts$N), ]

  if (ext == "png") {
    png(outfile, width = 2200, height = 1600, res = 250)
  } else {
    pdf(outfile, width = 8.8, height = 6.4)
  }

  par(mar = c(5, 9, 4, 2))
  bp <- barplot(
    counts$N,
    names.arg = counts$pattern,
    horiz = TRUE,
    las = 1,
    cex.names = 0.75,
    xlab = "Number of SNPs",
    main = paste0("Exclusive intersection counts: ", label, " ", mode)
  )

  text(
    x = counts$N,
    y = bp,
    labels = format(counts$N, big.mark = ",", scientific = FALSE),
    pos = 4,
    cex = 0.7
  )

  dev.off()
  message("Wrote: ", outfile)
}

for (label in c("1pct", "5pct")) {
  for (mode in c("projected", "strictHet")) {
    counts <- read_counts(label, mode)

    write.table(
      counts,
      file.path(OUTDIR, paste0("exclusive_counts_", label, "_", mode, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    plot_venn_like(counts, label, mode, "png")
    plot_venn_like(counts, label, mode, "pdf")
    plot_upset_bar(counts, label, mode, "png")
    plot_upset_bar(counts, label, mode, "pdf")
  }
}

message("All outputs written to: ", OUTDIR)
