suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
INDIR <- file.path(CAMP, "strict_HetDev_venn")
OUTDIR <- file.path(CAMP, "strict_HetDev_venn_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

ellipse_df <- function(cx, cy, rx, ry, angle_deg, set_name, n = 300) {
  t <- seq(0, 2*pi, length.out = n)
  a <- angle_deg * pi / 180
  x <- rx * cos(t)
  y <- ry * sin(t)
  data.table(
    x = cx + x * cos(a) - y * sin(a),
    y = cy + x * sin(a) + y * cos(a),
    set = set_name
  )
}

combo_name <- function(f, d, m, h) {
  labs <- c("FST", "DXY", "muLD", "HetDev")
  vals <- c(f, d, m, h)
  x <- labs[vals == 1]
  if (length(x) == 0) "none" else paste(x, collapse = "&")
}

make_counts <- function(dt) {
  dt[, FST_b := as_bool(FST)]
  dt[, DXY_b := as_bool(DXY)]
  dt[, muLD_b := as_bool(muLD)]
  dt[, Het_b := as_bool(HetDev_TX)]

  all_combos <- CJ(
    FST_i = 0:1,
    DXY_i = 0:1,
    muLD_i = 0:1,
    Het_i = 0:1
  )[!(FST_i == 0 & DXY_i == 0 & muLD_i == 0 & Het_i == 0)]

  dt[, key := paste0(as.integer(FST_b), as.integer(DXY_b), as.integer(muLD_b), as.integer(Het_b))]
  obs <- dt[, .N, by = key]

  all_combos[, key := paste0(FST_i, DXY_i, muLD_i, Het_i)]
  all_combos[, pattern := mapply(combo_name, FST_i, DXY_i, muLD_i, Het_i)]
  all_combos[, n_metrics := FST_i + DXY_i + muLD_i + Het_i]

  out <- merge(all_combos, obs, by = "key", all.x = TRUE)
  out[is.na(N), N := 0L]
  out[order(-n_metrics, pattern)]
}

plot_venn_like <- function(counts, label, mode) {
  ell <- rbindlist(list(
    ellipse_df(-0.75,  0.45, 1.25, 0.85,  30, "FST"),
    ellipse_df( 0.75,  0.45, 1.25, 0.85, -30, "DXY"),
    ellipse_df(-0.75, -0.45, 1.25, 0.85, -30, "muLD"),
    ellipse_df( 0.75, -0.45, 1.25, 0.85,  30, "HetDev_TX")
  ))

  pos <- data.table(
    key = c(
      "1000","0100","0010","0001",
      "1100","1010","1001","0110","0101","0011",
      "1110","1101","1011","0111","1111"
    ),
    x = c(
      -1.75, 1.75, -1.75, 1.75,
       0.00,-1.45, 0.00, 0.00, 1.45, 0.00,
      -0.45, 0.45,-0.45, 0.45, 0.00
    ),
    y = c(
       1.10, 1.10,-1.10,-1.10,
       1.25, 0.00, 0.35,-0.35, 0.00,-1.25,
       0.45, 0.45,-0.45,-0.45, 0.00
    )
  )

  lab <- merge(counts, pos, by = "key", all.x = TRUE)
  lab <- lab[N > 0]
  lab[, label_txt := paste0(pattern, "\n", N)]

  p <- ggplot() +
    geom_polygon(
      data = ell,
      aes(x = x, y = y, group = set, fill = set),
      alpha = 0.25,
      color = "grey30",
      linewidth = 0.6
    ) +
    geom_text(
      data = lab,
      aes(x = x, y = y, label = label_txt),
      size = 3
    ) +
    coord_equal(xlim = c(-2.6, 2.6), ylim = c(-1.9, 1.9)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0("Campbell SNP-level Venn-like diagram: ", label, " ", mode),
      subtitle = "Non-proportional schematic; labels show exclusive intersection counts",
      fill = "Metric"
    )

  safe <- paste0("Campbell_SNP_level_Venn_like_", label, "_", mode)

  ggsave(file.path(OUTDIR, paste0(safe, ".png")), p, width = 8, height = 6.5, dpi = 300)
  ggsave(file.path(OUTDIR, paste0(safe, ".pdf")), p, width = 8, height = 6.5)

  invisible(p)
}

plot_upset_bar <- function(counts, label, mode) {
  dt <- counts[N > 0][order(-n_metrics, -N)]
  dt[, pattern := factor(pattern, levels = rev(pattern))]

  p <- ggplot(dt, aes(x = pattern, y = N)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      title = paste0("Intersection counts: ", label, " ", mode),
      subtitle = "Exclusive SNP counts per pattern",
      x = "Intersection pattern",
      y = "Number of SNPs"
    )

  safe <- paste0("Campbell_SNP_level_UpSet_bar_", label, "_", mode)

  ggsave(file.path(OUTDIR, paste0(safe, ".png")), p, width = 8, height = 6.5, dpi = 300)
  ggsave(file.path(OUTDIR, paste0(safe, ".pdf")), p, width = 8, height = 6.5)

  invisible(p)
}

process_one <- function(label, mode) {
  file <- file.path(INDIR, paste0("membership_", label, "_", mode, ".tsv"))
  dt <- fread(file)

  counts <- make_counts(dt)

  fwrite(
    counts,
    file.path(OUTDIR, paste0("exclusive_counts_", label, "_", mode, ".tsv")),
    sep = "\t"
  )

  plot_venn_like(counts, label, mode)
  plot_upset_bar(counts, label, mode)

  cat("\n===", label, mode, "===\n")
  print(counts[N > 0][order(-n_metrics, -N)])
}

for (label in c("1pct", "5pct")) {
  for (mode in c("projected", "strictHet")) {
    process_one(label, mode)
  }
}

cat("\nOutputs written to:\n", OUTDIR, "\n")
