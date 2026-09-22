# ============================================================
# Regional plots for strict 1% FST+DXY+muLD Campbell-like regions
# Base R only
# ============================================================

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")

OUT <- file.path(CAMP, "Campbell_like_FST_DXY_muLD_candidate_region_plots_baseR")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

find_one <- function(pattern, root = CAMP) {
  cmd <- paste(
    "find",
    shQuote(root),
    "-type f | grep -E",
    shQuote(pattern),
    "| sort | head -1"
  )
  x <- system(cmd, intern = TRUE)
  if (length(x) == 0 || x[1] == "") {
    stop("Could not find file matching pattern: ", pattern)
  }
  x[1]
}

read_tab <- function(f) {
  message("Reading: ", f)
  if (grepl("\\.gz$", f)) {
    read.delim(gzfile(f), check.names = FALSE)
  } else {
    read.delim(f, check.names = FALSE)
  }
}

pick_col <- function(df, candidates, label) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(
      "No score column for ", label, ". Columns: ",
      paste(names(df), collapse = ", ")
    )
  }
  hit[1]
}

norm_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

# ------------------------------------------------------------
# Locate metric files
# ------------------------------------------------------------

FST_file <- find_one("FST_10kb.*windows.*\\.tsv\\.gz$")
DXY_file <- find_one("DXY_50SNP.*windows.*\\.tsv\\.gz$")
muLD_file <- find_one("muLD.*windows.*\\.tsv\\.gz$")

SNP1_file <- file.path(
  CAMP,
  "snp_level_campbell",
  "Campbell_SNP_membership_1pct_Campbell_threshold.tsv.gz"
)

if (!file.exists(SNP1_file)) {
  stop("Missing SNP membership file: ", SNP1_file)
}

# ------------------------------------------------------------
# Read metric files
# ------------------------------------------------------------

FST <- read_tab(FST_file)
DXY <- read_tab(DXY_file)
muLD <- read_tab(muLD_file)

fst_score <- pick_col(FST, c("FST_10kb_mean", "FST", "score"), "FST")
dxy_score <- pick_col(DXY, c("DXY_50SNP_mean", "DXY", "score"), "DXY")
muld_score <- pick_col(muLD, c("muLD_RAisd", "muLD", "LD", "mu", "score"), "muLD")

# Standardise core column names
standardise_metric <- function(df, score_col, metric_name) {
  chrom_col <- pick_col(df, c("CHROM", "chrom", "scaffold", "chr"), metric_name)
  start_col <- pick_col(df, c("window_start", "start", "lowPos", "BIN_START"), metric_name)
  end_col   <- pick_col(df, c("window_end", "end", "highPos", "BIN_END"), metric_name)

  out <- data.frame(
    CHROM = as.character(df[[chrom_col]]),
    start = as.numeric(df[[start_col]]),
    end = as.numeric(df[[end_col]]),
    score = as.numeric(df[[score_col]]),
    metric = metric_name,
    stringsAsFactors = FALSE
  )

  out <- out[is.finite(out$score) & !is.na(out$CHROM), ]
  out
}

FSTs <- standardise_metric(FST, fst_score, "FST")
DXYs <- standardise_metric(DXY, dxy_score, "DXY")
muLDs <- standardise_metric(muLD, muld_score, "muLD")

thresholds <- data.frame(
  metric = c("FST", "DXY", "muLD"),
  q95 = c(
    quantile(FSTs$score, 0.95, na.rm = TRUE),
    quantile(DXYs$score, 0.95, na.rm = TRUE),
    quantile(muLDs$score, 0.95, na.rm = TRUE)
  ),
  q99 = c(
    quantile(FSTs$score, 0.99, na.rm = TRUE),
    quantile(DXYs$score, 0.99, na.rm = TRUE),
    quantile(muLDs$score, 0.99, na.rm = TRUE)
  ),
  n_windows = c(nrow(FSTs), nrow(DXYs), nrow(muLDs)),
  score_col = c(fst_score, dxy_score, muld_score),
  stringsAsFactors = FALSE
)

write.table(
  thresholds,
  file.path(OUT, "FST_DXY_muLD_metric_thresholds_q95_q99.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Candidate regions
# ------------------------------------------------------------

regions <- data.frame(
  region_id = c(
    "CM044025.1_143590001_143600000",
    "CM044047.1_34830001_34850000",
    "CM044022.1_26030001_26040000",
    "CM044020.1_70690001_70700000",
    "CM044023.1_37750001_37760000"
  ),
  CHROM = c(
    "CM044025.1",
    "CM044047.1",
    "CM044022.1",
    "CM044020.1",
    "CM044023.1"
  ),
  start = c(143590001, 34830001, 26030001, 70690001, 37750001),
  end   = c(143600000, 34850000, 26040000, 70700000, 37760000),
  n_snps_expected = c(35, 30, 18, 5, 4),
  stringsAsFactors = FALSE
)

pad <- 100000
regions$plot_start <- pmax(1, regions$start - pad)
regions$plot_end <- regions$end + pad

write.table(
  regions,
  file.path(OUT, "FST_DXY_muLD_candidate_regions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Candidate SNPs from 1% membership
# ------------------------------------------------------------

SNP1 <- read_tab(SNP1_file)

SNP1$FST_b <- norm_bool(SNP1$FST)
SNP1$DXY_b <- norm_bool(SNP1$DXY)
SNP1$muLD_b <- norm_bool(SNP1$muLD)

cand_snps <- SNP1[
  SNP1$FST_b & SNP1$DXY_b & SNP1$muLD_b,
]

cand_snps$POS <- as.numeric(cand_snps$POS)

write.table(
  cand_snps,
  file.path(OUT, "candidate_SNPs_FST_DXY_muLD_1pct.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------

plot_region <- function(region_row) {

  chr <- region_row$CHROM
  ps <- region_row$plot_start
  pe <- region_row$plot_end
  rs <- region_row$start
  re <- region_row$end
  rid <- region_row$region_id

  metrics <- list(
    FST = FSTs,
    DXY = DXYs,
    muLD = muLDs
  )

  region_snps <- cand_snps[
    cand_snps$CHROM == chr &
      cand_snps$POS >= rs &
      cand_snps$POS <= re,
  ]

  snp_file <- file.path(OUT, paste0("candidate_SNPs_", rid, ".tsv"))
  write.table(
    region_snps,
    snp_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  png_file <- file.path(OUT, paste0(rid, "_regional_FST_DXY_muLD.png"))
  pdf_file <- file.path(OUT, paste0(rid, "_regional_FST_DXY_muLD.pdf"))

  for (device in c("png", "pdf")) {

    if (device == "png") {
      png(png_file, width = 1400, height = 1100, res = 150)
    } else {
      pdf(pdf_file, width = 10, height = 8)
    }

    par(mfrow = c(3, 1), mar = c(3.5, 4.5, 2.5, 1), oma = c(2, 0, 3, 0))

    for (metric_name in names(metrics)) {
      df <- metrics[[metric_name]]
      th <- thresholds[thresholds$metric == metric_name, ]

      dfr <- df[
        df$CHROM == chr &
          df$end >= ps &
          df$start <= pe,
      ]

      midpoint <- (dfr$start + dfr$end) / 2 / 1e6
      col <- ifelse(dfr$score >= th$q95, "darkblue", "lightblue")

      plot(
        midpoint,
        dfr$score,
        pch = 16,
        col = col,
        cex = 0.7,
        xlab = "",
        ylab = metric_name,
        main = paste0(metric_name, " windows"),
        xlim = c(ps, pe) / 1e6
      )

      abline(h = th$q95, lty = 2, lwd = 1.2)
      abline(h = th$q99, lty = 3, lwd = 1.2)

      abline(v = c(rs, re) / 1e6, lty = 3, col = "grey40")

      if (nrow(region_snps) > 0) {
        usr <- par("usr")
        rug_x <- region_snps$POS / 1e6
        segments(
          x0 = rug_x,
          x1 = rug_x,
          y0 = usr[3],
          y1 = usr[3] + 0.06 * (usr[4] - usr[3]),
          col = "black",
          lwd = 1
        )
      }

      legend(
        "topright",
        legend = c("Top 5% window", "Other window", "q95", "q99", "candidate SNPs"),
        col = c("darkblue", "lightblue", "black", "black", "black"),
        pch = c(16, 16, NA, NA, NA),
        lty = c(NA, NA, 2, 3, 1),
        bty = "n",
        cex = 0.75
      )
    }

    mtext(
      paste0(
        rid,
        " | strict 1% FST+DXY+muLD candidate region | n SNPs = ",
        nrow(region_snps)
      ),
      outer = TRUE,
      cex = 1.1,
      font = 2
    )

    mtext(
      paste0(chr, ":", format(rs, scientific = FALSE), "-", format(re, scientific = FALSE)),
      side = 1,
      outer = TRUE,
      cex = 0.9
    )

    dev.off()
  }

  data.frame(
    region_id = rid,
    CHROM = chr,
    start = rs,
    end = re,
    plot_start = ps,
    plot_end = pe,
    n_candidate_snps = nrow(region_snps),
    png = basename(png_file),
    pdf = basename(pdf_file),
    stringsAsFactors = FALSE
  )
}

plot_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(regions)), function(i) plot_region(regions[i, ]))
)

write.table(
  plot_summary,
  file.path(OUT, "FST_DXY_muLD_regional_plot_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nGenerated regional plots in:\n")
cat(OUT, "\n\n")
print(plot_summary)
