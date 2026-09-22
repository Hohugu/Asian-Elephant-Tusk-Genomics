BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
CORR <- file.path(CAMP, "venn_campbell_corrected")
CANDDIR <- file.path(CAMP, "strict_HetDev_venn")
OUTDIR <- file.path(CAMP, "Campbell_like_HetDev_candidate_region_plots_baseR")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (grepl("\\.gz$", path)) {
    read.delim(pipe(paste("zcat", shQuote(path))), stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

read_metric <- function(file, metric_name, score_candidates) {
  x <- read_tsv(file)
  score_col <- score_candidates[score_candidates %in% names(x)][1]

  if (is.na(score_col)) {
    stop(
      paste0(
        "No score column found for ", metric_name, " in ", file,
        "\nAvailable columns: ", paste(names(x), collapse = ", ")
      )
    )
  }

  needed <- c("CHROM", "window_start", "window_end")
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing columns in ", file, ": ",
        paste(missing, collapse = ", ")
      )
    )
  }

  out <- data.frame(
    CHROM = as.character(x$CHROM),
    window_start = as.numeric(x$window_start),
    window_end = as.numeric(x$window_end),
    score = as.numeric(x[[score_col]]),
    metric = metric_name,
    source_score_column = score_col,
    stringsAsFactors = FALSE
  )

  out$center <- (out$window_start + out$window_end) / 2
  out <- out[is.finite(out$score), ]
  out
}

message("Reading metrics...")

fst <- read_metric(
  file.path(CORR, "FST_10kb_windows.tsv.gz"),
  "FST",
  c("FST_10kb_mean", "FST", "FST_mean", "mean_FST", "weighted_FST")
)

dxy <- read_metric(
  file.path(CORR, "DXY_50SNP_step10_windows.tsv.gz"),
  "DXY",
  c("DXY_50SNP_mean", "DXY", "DXY_mean", "mean_DXY")
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

all_metrics <- rbind(fst, dxy, ld, het)

thr <- do.call(
  rbind,
  lapply(split(all_metrics, all_metrics$metric), function(z) {
    data.frame(
      metric = unique(z$metric),
      q95 = as.numeric(quantile(z$score, 0.95, na.rm = TRUE)),
      q99 = as.numeric(quantile(z$score, 0.99, na.rm = TRUE)),
      n_windows = nrow(z),
      source_score_column = unique(z$source_score_column)[1],
      stringsAsFactors = FALSE
    )
  })
)

write.table(
  thr,
  file.path(OUTDIR, "Campbell_like_metric_thresholds_q95_q99.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

make_regions <- function(cand, cluster_gap = 100000, flank = 100000) {
  cand$POS <- as.numeric(cand$POS)
  cand <- cand[order(cand$CHROM, cand$POS), ]

  region_id <- character(nrow(cand))
  for (chr in unique(cand$CHROM)) {
    ii <- which(cand$CHROM == chr)
    pos <- cand$POS[ii]
    ridx <- cumsum(c(TRUE, diff(pos) > cluster_gap))
    region_id[ii] <- paste0(chr, "_cluster_", ridx)
  }

  cand$region_id <- region_id

  regs <- do.call(
    rbind,
    lapply(split(cand, cand$region_id), function(z) {
      data.frame(
        region_id = unique(z$region_id),
        CHROM = unique(z$CHROM),
        region_start = min(z$POS),
        region_end = max(z$POS),
        n_snps = nrow(z),
        patterns = paste(sort(unique(z$pattern)), collapse = ";"),
        stringsAsFactors = FALSE
      )
    })
  )

  regs$plot_start <- pmax(1, regs$region_start - flank)
  regs$plot_end <- regs$region_end + flank
  regs$width_bp <- regs$region_end - regs$region_start + 1

  regs <- regs[order(regs$CHROM, regs$region_start), ]

  list(cand = cand, regions = regs)
}

plot_region <- function(cand, regions, mode_label, metrics_to_plot, suffix) {
  for (i in seq_len(nrow(regions))) {
    reg <- regions[i, ]

    snps_reg <- cand[cand$region_id == reg$region_id, ]
    snps_reg$pos_Mb <- snps_reg$POS / 1e6

    w <- all_metrics[
      all_metrics$metric %in% metrics_to_plot &
        all_metrics$CHROM == reg$CHROM &
        all_metrics$center >= reg$plot_start &
        all_metrics$center <= reg$plot_end,
    ]

    if (nrow(w) == 0) {
      warning("No metric windows found for ", reg$region_id, " ", mode_label)
      next
    }

    safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", paste(reg$region_id, mode_label, suffix, sep = "_"))

    for (ext in c("png", "pdf")) {
      outfile <- file.path(OUTDIR, paste0(safe_id, ".", ext))

      if (ext == "png") {
        png(outfile, width = 2500, height = ifelse(length(metrics_to_plot) == 3, 1800, 2200), res = 250)
      } else {
        pdf(outfile, width = 10, height = ifelse(length(metrics_to_plot) == 3, 7.2, 8.8))
      }

      oldpar <- par(no.readonly = TRUE)
      par(mfrow = c(length(metrics_to_plot), 1), mar = c(2.8, 4.5, 1.7, 1.5), oma = c(2.5, 0, 4.5, 0))

      for (metric_name in metrics_to_plot) {
        z <- w[w$metric == metric_name, ]
        z <- z[order(z$center), ]

        trow <- thr[thr$metric == metric_name, ]
        q95 <- trow$q95
        q99 <- trow$q99

        outlier <- z$score >= q95
        point_cols <- ifelse(outlier, "darkblue", "lightblue")

        plot(
          z$center / 1e6,
          z$score,
          pch = 16,
          cex = 0.75,
          col = point_cols,
          xlim = c(reg$plot_start, reg$plot_end) / 1e6,
          xlab = "",
          ylab = metric_name,
          main = metric_name
        )

        abline(h = q95, lty = 2, lwd = 1.2)
        abline(h = q99, lty = 3, lwd = 1.2)

        abline(v = snps_reg$POS / 1e6, col = adjustcolor("grey20", alpha.f = 0.35), lwd = 0.7)

        usr <- par("usr")
        rug(snps_reg$POS / 1e6, side = 1, col = "grey20", quiet = TRUE)
      }

      mtext(
        paste0(
          "Campbell-like magnified Manhattan: ",
          reg$CHROM, ":",
          format(reg$region_start, big.mark = ",", scientific = FALSE), "-",
          format(reg$region_end, big.mark = ",", scientific = FALSE),
          " (", mode_label, "; ", reg$n_snps, " SNPs)"
        ),
        outer = TRUE,
        side = 3,
        line = 2.5,
        cex = 1.1,
        font = 2
      )

      mtext(
        "Dashed = upper 5% quantile; dotted = upper 1% quantile; dark blue = upper 5% outlier windows; light blue = nonoutlier windows; vertical ticks = candidate SNPs",
        outer = TRUE,
        side = 3,
        line = 1.2,
        cex = 0.75
      )

      mtext(
        paste0(reg$CHROM, " position (Mb)"),
        outer = TRUE,
        side = 1,
        line = 0.8,
        cex = 0.95
      )

      par(oldpar)
      dev.off()
      message("Wrote: ", outfile)
    }
  }
}

process_mode <- function(cand_file, mode_label) {
  cand <- read.delim(cand_file, stringsAsFactors = FALSE, check.names = FALSE)

  reg_obj <- make_regions(cand)
  cand <- reg_obj$cand
  regs <- reg_obj$regions

  write.table(
    cand,
    file.path(OUTDIR, paste0("candidate_SNPs_", mode_label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    regs,
    file.path(OUTDIR, paste0("candidate_regions_", mode_label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  plot_region(
    cand,
    regs,
    mode_label,
    c("FST", "DXY", "muLD", "HetDev_TX"),
    "4panels_FST_DXY_muLD_Het"
  )

  plot_region(
    cand,
    regs,
    mode_label,
    c("FST", "muLD", "HetDev_TX"),
    "Campbell3panels_FST_muLD_Het"
  )

  message("Mode done: ", mode_label)
  message("SNPs: ", nrow(cand))
  message("Regions: ", nrow(regs))
}

projected_file <- file.path(CANDDIR, "SNPs_at_least3_including_HetDev_5pct_projected.tsv")
strict_file <- file.path(CANDDIR, "SNPs_at_least3_including_HetDev_5pct_strictHet.tsv")

process_mode(projected_file, "5pct_projected")
process_mode(strict_file, "5pct_strictHet")

message("All outputs written to: ", OUTDIR)
