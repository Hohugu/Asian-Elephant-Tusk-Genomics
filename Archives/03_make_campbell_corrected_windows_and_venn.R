camp <- Sys.getenv("CAMP")
if (camp == "") stop("CAMP is not set")

lib <- file.path(camp, "Rlibs")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", lib = lib, repos = "https://cloud.r-project.org")
}
library(data.table)

has_venn <- requireNamespace("VennDiagram", quietly = TRUE)
if (!has_venn) {
  try(
    install.packages("VennDiagram", lib = lib, repos = "https://cloud.r-project.org"),
    silent = TRUE
  )
  has_venn <- requireNamespace("VennDiagram", quietly = TRUE)
}

outdir <- file.path(camp, "venn_campbell_corrected")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fst_file <- file.path(camp, "metrics/fst/TT_vs_TX.weir.fst")
dxy_file <- file.path(camp, "metrics/dxy_het/DXY_HetDev_TX.full.tsv.gz")
ld_file  <- file.path(camp, "metrics/raisd/RAiSD_TX_muLD_full.tsv")

message("Input files:")
message(fst_file)
message(dxy_file)
message(ld_file)

# -----------------------------
# 1. FST: 10 kb non-overlapping windows
# -----------------------------
message("Reading FST")
fst <- fread(fst_file)
setnames(fst, "WEIR_AND_COCKERHAM_FST", "FST")
fst <- fst[!is.na(FST)]
fst[, POS := as.numeric(POS)]

message("FST: 10 kb windows")
fst[, window_start := floor((POS - 1) / 10000) * 10000 + 1]
fst[, window_end := window_start + 9999]

fst_w <- fst[, .(
  score = mean(FST, na.rm = TRUE),
  max_score = max(FST, na.rm = TRUE),
  n_snps = .N
), by = .(CHROM, window_start, window_end)]

fst_w[, metric := "FST"]
setnames(fst_w, "score", "FST_10kb_mean")

# -----------------------------
# 2. DXY and HetDev input
# -----------------------------
message("Reading DXY / HetDev")
dxy <- fread(
  cmd = paste("gzip -cd", shQuote(dxy_file)),
  select = c(
    "CHROM", "POS",
    "DXY",
    "HetDev_TX",
    "TX_p_alt",
    "TT_alt_sum",
    "Campbell_Het_candidate"
  )
)

dxy[, POS := as.numeric(POS)]
dxy <- dxy[!is.na(DXY)]

setorder(dxy, CHROM, POS)

# -----------------------------
# helper: SNP sliding windows
# -----------------------------
snp_sliding_windows <- function(dt, score_col, win_n, step_n, out_score_name) {
  dt <- dt[!is.na(get(score_col))]
  setorder(dt, CHROM, POS)

  res <- dt[, {
    if (.N < win_n) {
      data.table()
    } else {
      idx <- seq.int(1L, .N - win_n + 1L, by = step_n)
      vals <- get(score_col)

      means <- frollmean(vals, n = win_n, align = "left")[idx]

      data.table(
        window_start = POS[idx],
        window_end = POS[idx + win_n - 1L],
        score = means,
        n_snps = win_n
      )
    }
  }, by = CHROM]

  setnames(res, "score", out_score_name)
  res
}

# -----------------------------
# 3. DXY: 50-SNP sliding windows, step 10
# -----------------------------
message("DXY: 50-SNP sliding windows, step 10")
dxy_w <- snp_sliding_windows(
  dt = dxy,
  score_col = "DXY",
  win_n = 50,
  step_n = 10,
  out_score_name = "DXY_50SNP_mean"
)
dxy_w[, span_bp := window_end - window_start + 1]
dxy_w_span_filtered <- dxy_w[span_bp <= 10000]
dxy_w[, metric := "DXY"]
dxy_w_span_filtered[, metric := "DXY"]

# -----------------------------
# 4. Heterozygosity: TX 15-85%, absent in TT, 10-SNP step 2
# -----------------------------
message("HetDev: candidate SNP filter + 10-SNP sliding windows, step 2")

het_snps <- dxy[
  Campbell_Het_candidate == TRUE |
    (TX_p_alt >= 0.15 & TX_p_alt <= 0.85 & TT_alt_sum == 0)
]

het_snps <- het_snps[!is.na(HetDev_TX)]
setorder(het_snps, CHROM, POS)

het_w <- snp_sliding_windows(
  dt = het_snps,
  score_col = "HetDev_TX",
  win_n = 10,
  step_n = 2,
  out_score_name = "HetDev_TX_10SNP_mean"
)
het_w[, span_bp := window_end - window_start + 1]
het_w_span_filtered <- het_w[span_bp <= 10000]
het_w[, metric := "HetDev_TX"]
het_w_span_filtered[, metric := "HetDev_TX"]

# -----------------------------
# 5. LD: RAiSD muLD, default 50-SNP windows
# -----------------------------
message("Reading RAiSD muLD")
ld <- fread(ld_file, select = c("CHROM", "POS", "window_start", "window_end", "muLD"))
ld <- ld[!is.na(muLD)]
ld[, `:=`(
  POS = as.numeric(POS),
  window_start = as.numeric(window_start),
  window_end = as.numeric(window_end)
)]

ld_w <- ld[, .(
  CHROM,
  window_start,
  window_end,
  muLD_RAisd = muLD,
  n_snps = NA_integer_
)]
ld_w[, metric := "muLD"]

# -----------------------------
# Write method-specific window files
# -----------------------------
message("Writing method-specific window files")

fwrite(fst_w, file.path(outdir, "FST_10kb_windows.tsv.gz"), sep = "\t")
fwrite(dxy_w, file.path(outdir, "DXY_50SNP_step10_windows.tsv.gz"), sep = "\t")
fwrite(dxy_w_span_filtered, file.path(outdir, "DXY_50SNP_step10_windows_spanLE10kb.tsv.gz"), sep = "\t")
fwrite(het_w, file.path(outdir, "HetDev_TX_10SNP_step2_windows.tsv.gz"), sep = "\t")
fwrite(het_w_span_filtered, file.path(outdir, "HetDev_TX_10SNP_step2_windows_spanLE10kb.tsv.gz"), sep = "\t")
fwrite(ld_w,  file.path(outdir, "muLD_RAisd_50SNP_windows.tsv.gz"), sep = "\t")

window_counts <- data.table(
  metric = c("FST", "DXY", "muLD", "HetDev_TX"),
  method = c(
    "10kb non-overlapping windows",
    "50-SNP sliding windows, step 10 SNPs",
    "RAiSD muLD, default 50-SNP windows",
    "TX 15-85%, absent in TT; 10-SNP windows, step 2 SNPs"
  ),
  n_windows = c(nrow(fst_w), nrow(dxy_w), nrow(ld_w), nrow(het_w))
)
fwrite(window_counts, file.path(outdir, "campbell_corrected_window_counts.tsv"), sep = "\t")

# -----------------------------
# Normalize names for shared functions
# -----------------------------
fst_top_base <- fst_w[, .(
  metric = "FST",
  CHROM,
  window_start,
  window_end,
  score = FST_10kb_mean
)]

dxy_top_base <- dxy_w[, .(
  metric = "DXY",
  CHROM,
  window_start,
  window_end,
  score = DXY_50SNP_mean
)]

ld_top_base <- ld_w[, .(
  metric = "muLD",
  CHROM,
  window_start,
  window_end,
  score = muLD_RAisd
)]

het_top_base <- het_w[, .(
  metric = "HetDev_TX",
  CHROM,
  window_start,
  window_end,
  score = HetDev_TX_10SNP_mean
)]

bases <- list(
  FST = fst_top_base,
  DXY = dxy_top_base,
  muLD = ld_top_base,
  HetDev_TX = het_top_base
)

metrics <- names(bases)

# -----------------------------
# Convert selected windows to 10-kb genomic bins for overlap/Venn
# This is only for harmonizing different Campbell window definitions.
# -----------------------------
intervals_to_10kb_bins <- function(dt, bin_size = 10000) {
  dt <- copy(dt)
  dt[, window_start := as.numeric(window_start)]
  dt[, window_end := as.numeric(window_end)]
  dt[, bin_start_index := floor((window_start - 1) / bin_size)]
  dt[, bin_end_index := floor((window_end - 1) / bin_size)]
  dt[, row_id := .I]

  bins <- dt[, .(
    bin_index = seq.int(bin_start_index[1], bin_end_index[1])
  ), by = .(row_id, metric, CHROM)]

  bins[, bin_start := bin_index * bin_size + 1]
  bins[, bin_end := bin_start + bin_size - 1]
  bins[, bin_id := paste(CHROM, bin_start, bin_end, sep = ":")]

  unique(bins[, .(metric, CHROM, bin_start, bin_end, bin_id)])
}

make_threshold_outputs <- function(label, fraction) {
  message("Threshold: ", label, " fraction=", fraction)

  top_list <- list()

  for (m in metrics) {
    x <- copy(bases[[m]])
    x <- x[!is.na(score)]
    n_top <- max(1, ceiling(nrow(x) * fraction))

    x <- x[order(-score)]
    x <- x[seq_len(n_top)]
    x[, rank := seq_len(.N)]
    x[, threshold := label]

    top_list[[m]] <- x
  }

  top_windows <- rbindlist(top_list, use.names = TRUE)
  fwrite(
    top_windows,
    file.path(outdir, paste0("top_windows_", label, ".tsv.gz")),
    sep = "\t"
  )

  set_sizes_windows <- top_windows[, .N, by = .(threshold, metric)]
  setnames(set_sizes_windows, "N", "n_top_windows")
  fwrite(
    set_sizes_windows,
    file.path(outdir, paste0("set_sizes_top_windows_", label, ".tsv")),
    sep = "\t"
  )

  # Harmonized genomic bins for Venn
  bins <- intervals_to_10kb_bins(top_windows)
  fwrite(
    bins,
    file.path(outdir, paste0("top_windows_projected_to_10kb_bins_", label, ".tsv.gz")),
    sep = "\t"
  )

  sets <- lapply(metrics, function(m) bins[metric == m, unique(bin_id)])
  names(sets) <- metrics

  universe <- unique(unlist(sets, use.names = FALSE))
  mem <- data.table(bin_id = universe)

  for (m in metrics) {
    s <- sets[[m]]
    mem[, (m) := bin_id %in% s]
  }

  # recover genomic coordinates from bin_id
  parts <- tstrsplit(mem$bin_id, ":", fixed = TRUE)
  mem[, CHROM := parts[[1]]]
  mem[, bin_start := as.numeric(parts[[2]])]
  mem[, bin_end := as.numeric(parts[[3]])]
  setcolorder(mem, c("bin_id", "CHROM", "bin_start", "bin_end", metrics))

  fwrite(
    mem,
    file.path(outdir, paste0("venn_membership_10kb_bins_", label, ".tsv.gz")),
    sep = "\t"
  )

  pattern <- apply(as.matrix(mem[, ..metrics]), 1, function(v) {
    active <- metrics[as.logical(v)]
    if (length(active) == 0) "none" else paste(active, collapse = "&")
  })

  exact_counts <- data.table(pattern = pattern)[, .N, by = pattern][order(-N)]
  fwrite(
    exact_counts,
    file.path(outdir, paste0("venn_exact_intersection_counts_10kb_bins_", label, ".tsv")),
    sep = "\t"
  )

  set_sizes_bins <- data.table(
    threshold = label,
    metric = metrics,
    n_10kb_bins = sapply(metrics, function(m) sum(mem[[m]]))
  )
  fwrite(
    set_sizes_bins,
    file.path(outdir, paste0("venn_set_sizes_10kb_bins_", label, ".tsv")),
    sep = "\t"
  )

  pairwise <- CJ(metric1 = metrics, metric2 = metrics)
  pairwise <- pairwise[metric1 < metric2]
  pairwise[, n_10kb_bin_overlap := mapply(
    function(a, b) sum(mem[[a]] & mem[[b]]),
    metric1,
    metric2
  )]
  fwrite(
    pairwise,
    file.path(outdir, paste0("venn_pairwise_overlaps_10kb_bins_", label, ".tsv")),
    sep = "\t"
  )

  if (has_venn) {
    png_file <- file.path(outdir, paste0("Campbell_corrected_Venn_10kb_bins_", label, ".png"))

    VennDiagram::venn.diagram(
      x = sets,
      filename = png_file,
      imagetype = "png",
      height = 2600,
      width = 2600,
      resolution = 300,
      compression = "lzw",
      main = paste0("Campbell-style overlap, ", label),
      fill = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
      alpha = 0.35,
      cex = 0.65,
      cat.cex = 0.8,
      margin = 0.08
    )
  } else {
    message("VennDiagram unavailable; Venn counts written but no PNG.")
  }
}

# Campbell Fig. 2C uses 5% tails.
# Campbell Fig. 2D shows 5% and 1% thresholds.
thresholds <- data.table(
  label = c("5pct_Campbell", "1pct_Campbell_threshold"),
  fraction = c(0.05, 0.01)
)

for (i in seq_len(nrow(thresholds))) {
  make_threshold_outputs(thresholds$label[i], thresholds$fraction[i])
}

message("Done.")
message("Outputs in: ", outdir)

