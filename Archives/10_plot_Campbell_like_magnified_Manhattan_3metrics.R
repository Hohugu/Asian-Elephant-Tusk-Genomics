suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
CORR <- file.path(CAMP, "venn_campbell_corrected")
SNPDIR <- file.path(CAMP, "snp_level_campbell")
OUTDIR <- file.path(CAMP, "Campbell_like_regional_plots")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

membership_file <- file.path(SNPDIR, "Campbell_SNP_membership_5pct_Campbell.tsv.gz")

# =========================
# 1. SNPs in at least 3 metrics
# =========================

mem <- fread(
  membership_file,
  select = c(
    "snp_id", "CHROM", "POS", "REF", "ALT",
    "FST", "DXY", "muLD", "HetDev_TX",
    "n_analyses", "pattern",
    "DXY_value", "HetDev_TX_value", "TX_p_alt", "TT_alt_sum",
    "Campbell_Het_candidate"
  )
)

mem[, POS := as.numeric(POS)]
mem[, n_analyses := as.integer(n_analyses)]

cand <- mem[n_analyses >= 3]
setorder(cand, CHROM, POS)

fwrite(
  cand,
  file.path(OUTDIR, "Campbell_5pct_SNPs_at_least_3_metrics.tsv"),
  sep = "\t"
)

pattern_counts <- cand[, .N, by = pattern][order(-N)]
fwrite(
  pattern_counts,
  file.path(OUTDIR, "Campbell_5pct_SNPs_at_least_3_metrics_pattern_counts.tsv"),
  sep = "\t"
)

# =========================
# 2. Candidate regions / clusters
# =========================

cluster_gap <- 100000
plot_flank <- 100000

cand[, region_index := cumsum(c(TRUE, diff(POS) > cluster_gap)), by = CHROM]
cand[, region_id := paste0(CHROM, "_cluster_", region_index)]

regions <- cand[, .(
  region_start = min(POS),
  region_end = max(POS),
  n_snps = .N,
  patterns = paste(sort(unique(pattern)), collapse = ";")
), by = .(region_id, CHROM)]

regions[, plot_start := pmax(1, region_start - plot_flank)]
regions[, plot_end := region_end + plot_flank]
regions[, width_bp := region_end - region_start + 1]
setorder(regions, CHROM, region_start)

fwrite(
  regions,
  file.path(OUTDIR, "Campbell_5pct_regions_at_least_3_metrics.tsv"),
  sep = "\t"
)

# =========================
# 3. Load window-level metrics
# =========================

read_metric <- function(file, metric_name, score_candidates) {
  dt <- fread(file)

  score_col <- score_candidates[score_candidates %in% names(dt)][1]

  if (is.na(score_col)) {
    stop(
      paste0(
        "No score column found for ", metric_name, " in ", file,
        "\nAvailable columns: ", paste(names(dt), collapse = ", ")
      )
    )
  }

  needed <- c("CHROM", "window_start", "window_end")
  missing <- setdiff(needed, names(dt))
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing columns in ", file, ": ",
        paste(missing, collapse = ", "),
        "\nAvailable columns: ", paste(names(dt), collapse = ", ")
      )
    )
  }

  out <- dt[, .(
    CHROM = as.character(CHROM),
    window_start = as.numeric(window_start),
    window_end = as.numeric(window_end),
    score = as.numeric(get(score_col))
  )]

  out[, center := (window_start + window_end) / 2]
  out[, metric := metric_name]
  out[, source_score_column := score_col]

  out[is.finite(score)]
}

fst <- read_metric(
  file.path(CORR, "FST_10kb_windows.tsv.gz"),
  "FST",
  c("FST", "FST_mean", "mean_FST", "weighted_FST", "WEIR_AND_COCKERHAM_FST")
)

dxy <- read_metric(
  file.path(CORR, "DXY_50SNP_step10_windows.tsv.gz"),
  "DXY",
  c("DXY", "DXY_50SNP_mean", "DXY_mean", "mean_DXY")
)

ld <- read_metric(
  file.path(CORR, "muLD_RAisd_50SNP_windows.tsv.gz"),
  "muLD",
  c("muLD", "LD", "mu")
)

het <- read_metric(
  file.path(CORR, "HetDev_TX_10SNP_step2_windows.tsv.gz"),
  "HetDev_TX",
  c("HetDev_TX_10SNP_mean", "HetDev_TX", "HetDev", "score")
)

all_metrics <- rbindlist(list(fst, dxy, ld, het), use.names = TRUE)

# Upper 5% and 1% quantiles genome-wide, per metric
thr <- all_metrics[, .(
  q95 = as.numeric(quantile(score, 0.95, na.rm = TRUE)),
  q99 = as.numeric(quantile(score, 0.99, na.rm = TRUE)),
  n_windows = .N,
  source_score_column = unique(source_score_column)[1]
), by = metric]

fwrite(
  thr,
  file.path(OUTDIR, "Campbell_like_metric_thresholds_q95_q99.tsv"),
  sep = "\t"
)

# =========================
# 4. Plot each candidate region
# =========================

for (i in seq_len(nrow(regions))) {
  reg <- regions[i]

  snps_reg <- cand[region_id == reg$region_id]
  snps_reg[, pos_Mb := POS / 1e6]

  wdt <- all_metrics[
    CHROM == reg$CHROM &
      center >= reg$plot_start &
      center <= reg$plot_end
  ]

  if (nrow(wdt) == 0) {
    warning("No windows found for region: ", reg$region_id)
    next
  }

  wdt <- merge(wdt, thr[, .(metric, q95, q99)], by = "metric", all.x = TRUE)
  wdt[, outlier_5pct := score >= q95]
  wdt[, pos_Mb := center / 1e6]

  hlines <- unique(wdt[, .(metric, q95, q99)])

  title <- paste0(
    "Campbell-like magnified Manhattan: ",
    reg$CHROM, ":",
    format(reg$region_start, big.mark = ",", scientific = FALSE), "-",
    format(reg$region_end, big.mark = ",", scientific = FALSE)
  )

  subtitle <- paste0(
    reg$n_snps, " SNPs in >=3 metrics; ",
    "dashed = upper 5%, dotted = upper 1%; ",
    "vertical ticks = candidate SNPs"
  )

  p <- ggplot(wdt, aes(x = pos_Mb, y = score)) +
    geom_point(aes(color = outlier_5pct), size = 1.4, alpha = 0.85) +
    geom_hline(
      data = hlines,
      aes(yintercept = q95),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.35
    ) +
    geom_hline(
      data = hlines,
      aes(yintercept = q99),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_vline(
      data = snps_reg,
      aes(xintercept = pos_Mb),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.35
    ) +
    facet_grid(metric ~ ., scales = "free_y") +
    scale_color_manual(
      values = c("FALSE" = "#9ecae1", "TRUE" = "#084594"),
      labels = c("Nonoutlier window", "Upper 5% outlier window"),
      name = "Window status"
    ) +
    coord_cartesian(
      xlim = c(reg$plot_start, reg$plot_end) / 1e6
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text.y = element_text(angle = 0),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0(reg$CHROM, " position (Mb)"),
      y = "Window score"
    )

  safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", reg$region_id)

  ggsave(
    file.path(OUTDIR, paste0(safe_id, "_Campbell_like_magnified_Manhattan.pdf")),
    p,
    width = 9,
    height = 8
  )

  ggsave(
    file.path(OUTDIR, paste0(safe_id, "_Campbell_like_magnified_Manhattan.png")),
    p,
    width = 9,
    height = 8,
    dpi = 300
  )
}

cat("Done.\n")
cat("Output directory:", OUTDIR, "\n")
cat("SNPs >=3 metrics:", nrow(cand), "\n")
cat("Regions:", nrow(regions), "\n")
