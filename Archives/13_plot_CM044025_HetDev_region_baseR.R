B <- "/scratch/project_2000886/Hoedric/GWAS_2025"
C <- file.path(B, "Genetics_Analysis/Campbell_exact_reproduction")
D <- file.path(C, "venn_campbell_corrected")
S <- file.path(C, "strict_HetDev_venn")
O <- file.path(C, "Campbell_like_HetDev_candidate_region_plots_baseR")
dir.create(O, recursive = TRUE, showWarnings = FALSE)

rd <- function(f) {
  if (!file.exists(f)) stop("Missing file: ", f)
  if (grepl("\\.gz$", f)) {
    read.delim(pipe(paste("zcat", shQuote(f))), check.names = FALSE)
  } else {
    read.delim(f, check.names = FALSE)
  }
}

pick <- function(x) {
  for (f in x) if (file.exists(f)) return(f)
  stop("No file found among: ", paste(x, collapse = ", "))
}

metric <- function(f, m, cols) {
  x <- rd(f)
  cc <- cols[cols %in% names(x)][1]
  if (is.na(cc)) stop("No score column for ", m, ". Columns: ", paste(names(x), collapse = ", "))
  y <- data.frame(
    CHROM = as.character(x$CHROM),
    start = as.numeric(x$window_start),
    end = as.numeric(x$window_end),
    score = as.numeric(x[[cc]]),
    metric = m,
    score_col = cc
  )
  y$center <- (y$start + y$end) / 2
  y[is.finite(y$score), ]
}

message("Reading metrics")

F <- metric(file.path(D, "FST_10kb_windows.tsv.gz"), "FST",
            c("FST_10kb_mean", "FST", "FST_mean", "mean_FST", "score"))

X <- metric(file.path(D, "DXY_50SNP_step10_windows.tsv.gz"), "DXY",
            c("DXY_50SNP_mean", "DXY", "DXY_mean", "mean_DXY", "score"))

L <- metric(file.path(D, "muLD_RAisd_50SNP_windows.tsv.gz"), "muLD",
            c("muLD_RAisd", "muLD", "LD", "mu", "score"))

H <- metric(pick(c(
              file.path(D, "HetDev_TX_10SNP_step2_windows_spanLE10kb.tsv.gz"),
              file.path(D, "HetDev_TX_10SNP_step2_windows.tsv.gz")
            )), "HetDev_TX",
            c("HetDev_TX_10SNP_mean", "HetDev_TX", "HetDev", "score"))

M <- rbind(F, X, L, H)

TH <- do.call(rbind, lapply(split(M, M$metric), function(z) {
  data.frame(
    metric = unique(z$metric),
    q95 = as.numeric(quantile(z$score, 0.95, na.rm = TRUE)),
    q99 = as.numeric(quantile(z$score, 0.99, na.rm = TRUE)),
    n_windows = nrow(z),
    score_col = unique(z$score_col)[1]
  )
}))

write.table(TH, file.path(O, "Campbell_like_metric_thresholds_q95_q99.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

P <- rd(file.path(S, "SNPs_at_least3_including_HetDev_5pct_projected.tsv"))
Q <- rd(file.path(S, "SNPs_at_least3_including_HetDev_5pct_strictHet.tsv"))

write.table(P, file.path(O, "candidate_SNPs_5pct_projected.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(Q, file.path(O, "candidate_SNPs_5pct_strictHet.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

chr <- "CM044025.1"
a <- 22575792
b <- 22576805
ps <- a - 100000
pe <- b + 100000

REG <- data.frame(
  mode = c("5pct_projected", "5pct_strictHet"),
  CHROM = chr,
  region_start = a,
  region_end = b,
  n_snps = c(nrow(P), nrow(Q)),
  plot_start = ps,
  plot_end = pe,
  width_bp = b - a + 1
)

write.table(REG, file.path(O, "candidate_regions_CM044025_HetDev.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plotset <- function(snps, lab, mets, suffix) {
  snps$POS <- as.numeric(snps$POS)
  W <- M[M$metric %in% mets & M$CHROM == chr & M$center >= ps & M$center <= pe, ]
  if (nrow(W) == 0) stop("No metric windows in plotting region")

  id <- paste0("CM044025.1_22575792_22576805_", lab, "_", suffix)

  for (ext in c("png", "pdf")) {
    out <- file.path(O, paste0(id, ".", ext))

    if (ext == "png") {
      png(out, width = 2500, height = ifelse(length(mets) == 3, 1800, 2200), res = 250)
    } else {
      pdf(out, width = 10, height = ifelse(length(mets) == 3, 7.2, 8.8))
    }

    old <- par(no.readonly = TRUE)
    par(mfrow = c(length(mets), 1),
        mar = c(2.8, 4.7, 1.6, 1.2),
        oma = c(2.4, 0, 4.6, 0))

    for (m in mets) {
      z <- W[W$metric == m, ]
      z <- z[order(z$center), ]
      t <- TH[TH$metric == m, ]

      plot(z$center / 1e6, z$score,
           pch = 16,
           cex = 0.75,
           col = ifelse(z$score >= t$q95, "darkblue", "lightblue"),
           xlim = c(ps, pe) / 1e6,
           xlab = "",
           ylab = m,
           main = m)

      abline(h = t$q95, lty = 2, lwd = 1.2)
      abline(h = t$q99, lty = 3, lwd = 1.2)
      abline(v = snps$POS / 1e6,
             col = adjustcolor("grey20", alpha.f = 0.35),
             lwd = 0.7)
      rug(snps$POS / 1e6, side = 1, col = "grey20", quiet = TRUE)
    }

    mtext(paste0("Campbell-like magnified Manhattan: ", chr, ":",
                 format(a, big.mark = ",", scientific = FALSE), "-",
                 format(b, big.mark = ",", scientific = FALSE),
                 " (", lab, "; ", nrow(snps), " SNPs)"),
          outer = TRUE, side = 3, line = 2.6, cex = 1.1, font = 2)

    mtext("Dashed = upper 5%; dotted = upper 1%; dark blue = upper 5% outlier windows; light blue = nonoutlier windows; vertical ticks = candidate SNPs",
          outer = TRUE, side = 3, line = 1.3, cex = 0.75)

    mtext(paste0(chr, " position (Mb)"),
          outer = TRUE, side = 1, line = 0.7, cex = 0.95)

    par(old)
    dev.off()
    message("Wrote: ", out)
  }
}

plotset(P, "5pct_projected",
        c("FST", "DXY", "muLD", "HetDev_TX"),
        "4panels_FST_DXY_muLD_Het")

plotset(P, "5pct_projected",
        c("FST", "muLD", "HetDev_TX"),
        "Campbell3panels_FST_muLD_Het")

plotset(Q, "5pct_strictHet",
        c("FST", "DXY", "muLD", "HetDev_TX"),
        "4panels_FST_DXY_muLD_Het")

plotset(Q, "5pct_strictHet",
        c("FST", "muLD", "HetDev_TX"),
        "Campbell3panels_FST_muLD_Het")

message("All outputs written to: ", O)
