#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
options(scipen = 999)

# ============================================================
# 25bis_plot_genomewide_scans_with_GEMMA_Campbell.R
#
# Re-plot genome-wide population-genomic scans with:
#   - GEMMA all-sample Bonferroni SNPs
#   - GEMMA male-only Bonferroni SNPs
#   - Campbell candidate genes
#
# No recombination layer is used here.
# ============================================================

BASE <- Sys.getenv("BASE", unset = "/scratch/project_2000886/Hoedric/GWAS_2025")

# ---- Input paths ----
FST_50KB <- file.path(BASE, "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_windows.tsv")
HO_50KB  <- file.path(BASE, "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv")
PI_50KB  <- file.path(BASE, "Genetics_Analysis/Diversity/TT_TX_pi_50kb_comparison.tsv")

TAJ_CANDIDATES <- c(
  file.path(BASE, "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_comparison.tsv"),
  file.path(BASE, "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_windows.tsv"),
  file.path(BASE, "Genetics_Analysis/Diversity/TT_TX_TajimaD_comparison.tsv")
)
TAJ_50KB <- TAJ_CANDIDATES[file.exists(TAJ_CANDIDATES)][1]
if (is.na(TAJ_50KB)) {
  found <- list.files(file.path(BASE, "Genetics_Analysis/Diversity"),
                      pattern = "TajimaD.*50kb.*(comparison|window).*\\.(tsv|txt)$",
                      full.names = TRUE, ignore.case = TRUE)
  TAJ_50KB <- found[1]
}

GEMMA_ALL <- file.path(BASE, "GWAS/result/GEMMA_PC3/output/GEMMA_PC3_Bonferroni_hits.tsv")
if (!file.exists(GEMMA_ALL)) {
  GEMMA_ALL <- file.path(BASE, "GWAS/result/GEMMA_PC3/output/GEMMA_PC3_Bonferroni_hits.txt")
}
GEMMA_MALE <- file.path(BASE, "GWAS/tables_GEMMA_male/GEMMA_male_Bonferroni_hits.tsv")
if (!file.exists(GEMMA_MALE)) {
  GEMMA_MALE <- file.path(BASE, "GWAS/result/GEMMA_male/output/GEMMA_male_PC3_Bonferroni_hits.txt")
}

GFF <- file.path(BASE, "../Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff")
if (!file.exists(GFF)) {
  GFF <- "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"
}
MAP <- file.path(BASE, "GWAS/tables_GEMMA_male/CM_to_NC.map")

OUTDIR <- file.path(BASE, "Genetics_Analysis/Diversity/figures/GEMMA_Campbell_annotated")
TABDIR <- file.path(BASE, "Genetics_Analysis/Diversity/tables")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Plot settings ----
LABEL_ALL_GEMMA <- Sys.getenv("LABEL_ALL_GEMMA", unset = "0") == "1"
N_LABEL_ALL_TOP <- as.integer(Sys.getenv("N_LABEL_ALL_TOP", unset = "25"))
MAIN_WIDTH  <- 3400
MAIN_HEIGHT <- 1250
PNG_RES     <- 170

COL_BACKGROUND <- "grey75"
COL_TOP        <- "orange"
COL_EXTREME    <- "blue"
COL_GEMMA_ALL  <- "red3"
COL_GEMMA_MALE <- "purple4"
COL_CAMPBELL   <- "darkgreen"
COL_LINE       <- "grey88"

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
stop_missing <- function(path, label) {
  if (is.na(path) || !file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

stop_missing(FST_50KB, "FST 50 kb table")
stop_missing(HO_50KB,  "Ho 50 kb table")
stop_missing(PI_50KB,  "pi 50 kb table")
stop_missing(TAJ_50KB, "TajimaD 50 kb table")
stop_missing(GEMMA_ALL, "GEMMA all-sample Bonferroni table")
stop_missing(GEMMA_MALE, "GEMMA male-only Bonferroni table")

# ============================================================
# Helper functions
# ============================================================

cm_chr_num <- function(x) {
  suppressWarnings(as.integer(sub("^CM0440", "", sub("\\.1$", "", x))))
}

keep_cm20_48 <- function(dt, chrom_col = "CHROM") {
  dt <- copy(dt)
  dt[, chr_num := cm_chr_num(get(chrom_col))]
  dt[!is.na(chr_num) & chr_num >= 20 & chr_num <= 48]
}

first_existing_col <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) NA_character_ else hit
}

standardize_windows <- function(dt, chrom_col, start_col, end_col, value_col,
                                snp_col = NULL, min_snps = NULL) {
  d <- copy(dt)
  setnames(d, chrom_col, "CHROM")
  setnames(d, start_col, "START")
  setnames(d, end_col, "END")
  setnames(d, value_col, "VALUE")
  d[, CHROM := as.character(CHROM)]
  d[, START := as.numeric(START)]
  d[, END := as.numeric(END)]
  d[, VALUE := as.numeric(VALUE)]
  if (!is.null(snp_col) && snp_col %in% names(d) && !is.null(min_snps)) {
    d <- d[as.numeric(get(snp_col)) >= min_snps]
  }
  d <- keep_cm20_48(d, "CHROM")
  d <- d[!is.na(VALUE) & !is.na(START) & !is.na(END)]
  setorder(d, chr_num, START)
  d
}

find_tajima_value_col <- function(dt) {
  candidates <- c("delta_TajimaD", "delta_TajimaD_TX_minus_TT", "TajimaD_TX_minus_TT",
                  "delta_D", "Delta_TajimaD", "tajimaD_delta")
  c1 <- first_existing_col(dt, candidates)
  if (!is.na(c1)) return(c1)
  numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("BIN_START", "BIN_END", "START", "END", "N_VARIANTS", "N_SNPs"))
  if (length(numeric_cols) == 0) stop("Could not detect TajimaD delta column", call. = FALSE)
  numeric_cols[1]
}

make_layout <- function(...) {
  allw <- rbindlist(list(...), fill = TRUE)
  allw <- allw[!is.na(chr_num) & !is.na(END)]
  chr_sizes <- allw[, .(chr_len = as.numeric(max(END, na.rm = TRUE))), by = .(CHROM, chr_num)]
  setorder(chr_sizes, chr_num)
  chr_sizes[, offset := c(0, cumsum(as.numeric(chr_len))[1:(.N - 1)])]
  chr_sizes[, center := offset + chr_len / 2]
  chr_sizes[, label := as.character(chr_num)]
  chr_sizes
}

add_xpos_windows <- function(d, layout) {
  d <- merge(d, layout[, .(CHROM, chr_num, offset)], by = c("CHROM", "chr_num"), all.x = TRUE)
  d[, xpos := START + offset]
  setorder(d, chr_num, START)
  d
}

parse_snp_id_vector <- function(x) {
  # Expected form: CM044048.1:14979125:G:T
  s <- tstrsplit(as.character(x), ":", fixed = TRUE, keep = 1:2)
  data.table(
    SNP_ID = as.character(x),
    CHROM = s[[1]],
    POS = suppressWarnings(as.numeric(s[[2]]))
  )
}

detect_snp_col <- function(dt) {
  candidates <- c("SNP", "snp", "rs", "ID", "MarkerName", "variant", "Variant", "SNP_ID")
  for (cc in candidates) {
    if (cc %in% names(dt)) {
      vals <- as.character(dt[[cc]])
      if (sum(grepl("^CM0440[0-9]+\\.1:[0-9]+", vals), na.rm = TRUE) > 0) return(cc)
    }
  }
  for (cc in names(dt)) {
    vals <- as.character(dt[[cc]])
    if (sum(grepl("^CM0440[0-9]+\\.1:[0-9]+", vals), na.rm = TRUE) > 0) return(cc)
  }
  stop("Could not detect SNP column in GEMMA table", call. = FALSE)
}

detect_p_col <- function(dt) {
  candidates <- c("p_wald", "P", "p", "p_lrt", "p_score", "pvalue", "p_value", "P_VALUE", "p_wald_all")
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

read_gemma_hits <- function(path, model_name) {
  msg("Reading GEMMA hits:", path)
  dt <- fread(path)
  snp_col <- detect_snp_col(dt)
  p_col <- detect_p_col(dt)
  snps <- parse_snp_id_vector(dt[[snp_col]])
  snps[, model := model_name]
  if (!is.na(p_col)) snps[, P := suppressWarnings(as.numeric(dt[[p_col]]))] else snps[, P := NA_real_]
  snps <- keep_cm20_48(snps, "CHROM")
  snps <- snps[!is.na(POS)]
  snps[, short_label := paste0(ifelse(model == "GEMMA_all_PC3", "All", "Male"),
                                "_", chr_num, ":", sprintf("%.2f", POS / 1e6))]
  unique(snps, by = c("SNP_ID", "model"))
}

read_cm_nc_map <- function(path) {
  if (!file.exists(path)) return(NULL)
  m <- tryCatch(fread(path, header = FALSE, fill = TRUE), error = function(e) NULL)
  if (is.null(m) || ncol(m) < 2) return(NULL)
  cm_col <- which(vapply(m, function(z) any(grepl("^CM", as.character(z))), logical(1)))[1]
  nc_col <- which(vapply(m, function(z) any(grepl("^NC", as.character(z))), logical(1)))[1]
  if (is.na(cm_col) || is.na(nc_col)) return(NULL)
  unique(data.table(CHROM = as.character(m[[cm_col]]), NC_chr = as.character(m[[nc_col]]))[
    grepl("^CM", CHROM) & grepl("^NC", NC_chr)
  ])
}

extract_attr <- function(attr, key) {
  # Extracts key=value until semicolon
  pat <- paste0("(^|;)", key, "=([^;]+)")
  out <- sub(pat, "\\2", attr)
  out[!grepl(pat, attr)] <- NA_character_
  out
}

read_campbell_from_gff <- function(gff, map_path) {
  symbols <- c("AMELX", "AMBN", "AMTN", "ENAM", "ODAM", "MEP1A", "MEP1B", "PLA2G7")
  if (!file.exists(gff)) return(data.table())

  msg("Reading candidate genes from GFF:", gff)
  cmd <- paste("awk '$3==\"gene\"'", shQuote(gff))
  g <- tryCatch(fread(cmd = cmd, sep = "\t", header = FALSE, fill = TRUE, quote = ""),
                error = function(e) data.table())
  if (nrow(g) == 0 || ncol(g) < 9) return(data.table())
  setnames(g, paste0("V", seq_len(ncol(g))))
  g <- g[, .(seqid = as.character(V1), feature = V3, START = as.numeric(V4),
             END = as.numeric(V5), strand = V7, attr = as.character(V9))]

  g[, gene_id := fifelse(!is.na(extract_attr(attr, "ID")), extract_attr(attr, "ID"), extract_attr(attr, "Dbxref"))]
  g[, Name := extract_attr(attr, "Name")]
  g[, gene := extract_attr(attr, "gene")]

  hits <- rbindlist(lapply(symbols, function(sym) {
    x <- g[grepl(paste0("(^|[^A-Za-z0-9])", sym, "([^A-Za-z0-9]|$)"), attr)]
    if (nrow(x) == 0) return(data.table())
    x[, Campbell_gene := sym]
    x
  }), fill = TRUE)

  if (nrow(hits) == 0) return(data.table())

  hits[, CHROM := seqid]
  mp <- read_cm_nc_map(map_path)
  if (!is.null(mp)) {
    hits <- merge(hits, mp, by.x = "seqid", by.y = "NC_chr", all.x = TRUE, suffixes = c(".gff", ".map"))
    if ("CHROM.map" %in% names(hits)) {
      hits[, CHROM := fifelse(!is.na(`CHROM.map`), `CHROM.map`, `CHROM.gff`)]
      hits[, c("CHROM.gff", "CHROM.map") := NULL]
    }
  }

  hits <- hits[grepl("^CM0440", CHROM)]
  hits[, chr_num := cm_chr_num(CHROM)]
  hits <- hits[chr_num >= 20 & chr_num <= 48]
  hits[, POS := (START + END) / 2]
  hits[, marker_type := "Campbell_gene"]
  hits[, SNP_ID := NA_character_]
  hits[, P := NA_real_]
  hits[, model := "Campbell_candidate_gene"]
  hits[, short_label := Campbell_gene]
  hits[, gene_label := Campbell_gene]

  # Keep one row per Campbell symbol, preferring exact gene field when possible and longest interval otherwise.
  hits[, exact_gene := fifelse(!is.na(gene) & gene == Campbell_gene, 1L, 0L)]
  hits[, len := END - START]
  setorder(hits, Campbell_gene, -exact_gene, -len)
  hits <- hits[, .SD[1], by = Campbell_gene]
  hits[, .(SNP_ID, model, marker_type, Campbell_gene, gene_label, CHROM, START, END, POS, P, short_label, chr_num)]
}

fallback_campbell <- function() {
  data.table(
    Campbell_gene = c("ODAM", "MEP1A", "AMTN", "ENAM", "AMELX", "AMBN"),
    CHROM = c("CM044024.1", "CM044020.1", "CM044024.1", "CM044024.1", "CM044047.1", "CM044024.1"),
    START = c(87832910, 111679947, 87416785, 87316884, 168725165, 87354290),
    END   = c(87841335, 111709525, 87432258, 87332008, 168729418, 87366652)
  )[, `:=`(
    POS = (START + END) / 2,
    SNP_ID = NA_character_,
    P = NA_real_,
    model = "Campbell_candidate_gene",
    marker_type = "Campbell_gene",
    gene_label = Campbell_gene,
    short_label = Campbell_gene,
    chr_num = cm_chr_num(CHROM)
  )]
}

attach_layout_to_markers <- function(markers, layout) {
  m <- copy(markers)
  m[, chr_num := cm_chr_num(CHROM)]
  m <- merge(m, layout[, .(CHROM, chr_num, offset)], by = c("CHROM", "chr_num"), all.x = TRUE)
  m[, xpos := POS + offset]
  m[!is.na(xpos)]
}

value_at_markers <- function(windows, markers) {
  if (nrow(markers) == 0) return(markers[, metric_value := numeric()])
  out <- copy(markers)
  out[, metric_value := NA_real_]
  for (i in seq_len(nrow(out))) {
    chr <- out$CHROM[i]
    pos <- out$POS[i]
    w <- windows[CHROM == chr & START <= pos & END >= pos]
    if (nrow(w) == 0) {
      w <- windows[CHROM == chr]
      if (nrow(w) > 0) {
        w[, center_tmp := (START + END) / 2]
        w <- w[which.min(abs(center_tmp - pos))]
      }
    }
    if (nrow(w) > 0) out$metric_value[i] <- w$VALUE[1]
  }
  out
}

label_gemma_subset <- function(gemma_all, gemma_male) {
  # Label all male-only markers and only top all-sample markers by p-value by default.
  male_lab <- copy(gemma_male)
  if (LABEL_ALL_GEMMA) {
    all_lab <- copy(gemma_all)
  } else {
    all_lab <- copy(gemma_all[order(P)])
    all_lab <- all_lab[seq_len(min(N_LABEL_ALL_TOP, nrow(all_lab)))]
  }
  unique(rbindlist(list(all_lab, male_lab), fill = TRUE), by = c("SNP_ID", "model"))
}

draw_text_staggered <- function(x, y, labels, col, y_range, cex = 0.50, srt = 50) {
  if (length(x) == 0) return(invisible(NULL))
  ok <- !is.na(x) & !is.na(y) & !is.na(labels)
  x <- x[ok]; y <- y[ok]; labels <- labels[ok]
  if (length(x) == 0) return(invisible(NULL))
  ord <- order(x)
  x <- x[ord]; y <- y[ord]; labels <- labels[ord]
  step <- 0.045 * diff(y_range)
  ytext <- y + step * (1 + ((seq_along(y) - 1) %% 4))
  ytext <- pmin(ytext, y_range[2] - 0.02 * diff(y_range))
  segments(x, y, x, ytext, col = col, lwd = 0.4)
  text(x, ytext, labels = labels, col = col, cex = cex, srt = srt, adj = c(0, 0.5), xpd = NA)
}

draw_campbell_lines <- function(campbell, y_range) {
  if (nrow(campbell) == 0) return(invisible(NULL))
  campbell <- campbell[!is.na(xpos)]
  if (nrow(campbell) == 0) return(invisible(NULL))
  abline(v = campbell$xpos, col = COL_CAMPBELL, lty = 3, lwd = 1.2)
  yr <- diff(y_range)
  top <- y_range[2] - 0.015 * yr
  campbell <- campbell[order(xpos)]
  for (i in seq_len(nrow(campbell))) {
    yy <- top - ((i - 1) %% 5) * 0.045 * yr
    text(campbell$xpos[i], yy, labels = campbell$gene_label[i], srt = 90,
         col = COL_CAMPBELL, cex = 0.62, adj = c(1, 0.5), xpd = NA)
  }
}

plot_scan <- function(windows, layout, gemma_all, gemma_male, campbell,
                      outfile_prefix, main, ylab, mode = c("high", "low", "two_sided")) {
  mode <- match.arg(mode)

  d <- copy(windows)
  y <- d$VALUE
  q01  <- as.numeric(quantile(y, 0.01,  na.rm = TRUE))
  q001 <- as.numeric(quantile(y, 0.001, na.rm = TRUE))
  q99  <- as.numeric(quantile(y, 0.99,  na.rm = TRUE))
  q999 <- as.numeric(quantile(y, 0.999, na.rm = TRUE))

  d[, category := "background"]
  if (mode == "high") {
    d[VALUE >= q99,  category := "top_1"]
    d[VALUE >= q999, category := "top_0.1"]
    threshold_lines <- data.table(y = c(q99, q999), col = c(COL_TOP, COL_EXTREME), lty = 2)
    legend_metric <- c("Background", "top 1%", "top 0.1%")
    legend_metric_col <- c(COL_BACKGROUND, COL_TOP, COL_EXTREME)
  } else if (mode == "low") {
    d[VALUE <= q01,  category := "bottom_1"]
    d[VALUE <= q001, category := "bottom_0.1"]
    threshold_lines <- data.table(y = c(q01, q001), col = c(COL_TOP, COL_EXTREME), lty = 2)
    legend_metric <- c("Background", "bottom 1%", "bottom 0.1%")
    legend_metric_col <- c(COL_BACKGROUND, COL_TOP, COL_EXTREME)
  } else {
    d[VALUE <= q01, category := "bottom_1"]
    d[VALUE >= q99, category := "top_1"]
    threshold_lines <- data.table(y = c(q01, q99), col = c(COL_EXTREME, COL_TOP), lty = 2)
    legend_metric <- c("Background", "bottom 1%", "top 1%")
    legend_metric_col <- c(COL_BACKGROUND, COL_EXTREME, COL_TOP)
  }

  d[, point_col := COL_BACKGROUND]
  d[category %in% c("top_1", "bottom_1"), point_col := COL_TOP]
  d[category %in% c("top_0.1", "bottom_0.1"), point_col := COL_EXTREME]
  d[, point_cex := fifelse(category %in% c("top_0.1", "bottom_0.1"), 0.55,
                          fifelse(category %in% c("top_1", "bottom_1"), 0.45, 0.25))]

  ga <- value_at_markers(d, gemma_all)
  gm <- value_at_markers(d, gemma_male)
  label_gemma <- value_at_markers(d, label_gemma_subset(gemma_all, gemma_male))

  y_range <- range(c(d$VALUE, ga$metric_value, gm$metric_value), na.rm = TRUE)
  if (diff(y_range) == 0) y_range <- y_range + c(-1, 1)
  y_pad <- 0.18 * diff(y_range)
  y_range <- y_range + c(-0.08 * diff(y_range), y_pad)

  x_range <- range(c(layout$offset, layout$offset + layout$chr_len), na.rm = TRUE)
  chr_boundaries <- layout$offset[-1]

  png_file <- file.path(OUTDIR, paste0(outfile_prefix, ".png"))
  pdf_file <- file.path(OUTDIR, paste0(outfile_prefix, ".pdf"))

  plot_fun <- function() {
    par(mar = c(6.2, 5.2, 4.2, 2.0), cex.main = 1.15, cex.lab = 1.05, cex.axis = 0.75)
    plot(d$xpos, d$VALUE,
         pch = 20, cex = d$point_cex, col = d$point_col,
         xaxt = "n", xlab = "Chromosome", ylab = ylab,
         main = main, xlim = x_range, ylim = y_range)
    abline(v = chr_boundaries, col = COL_LINE, lwd = 0.45)
    axis(1, at = layout$center, labels = layout$label, las = 2, cex.axis = 0.72)
    for (i in seq_len(nrow(threshold_lines))) {
      abline(h = threshold_lines$y[i], col = threshold_lines$col[i], lty = threshold_lines$lty[i], lwd = 1.4)
    }
    if (min(y_range) < 0 && max(y_range) > 0) abline(h = 0, col = "grey40", lty = 2, lwd = 0.8)

    # Campbell genes as vertical lines
    draw_campbell_lines(campbell, y_range)

    # GEMMA overlays
    if (nrow(ga) > 0) points(ga$xpos, ga$metric_value, pch = 24, bg = COL_GEMMA_ALL,
                              col = COL_GEMMA_ALL, cex = 0.78, lwd = 0.8)
    if (nrow(gm) > 0) points(gm$xpos, gm$metric_value, pch = 23, bg = COL_GEMMA_MALE,
                              col = COL_GEMMA_MALE, cex = 0.86, lwd = 0.9)

    draw_text_staggered(label_gemma$xpos, label_gemma$metric_value, label_gemma$short_label,
                        col = ifelse(label_gemma$model == "GEMMA_male_PC3", COL_GEMMA_MALE, COL_GEMMA_ALL),
                        y_range = y_range, cex = 0.46, srt = 45)

    legend("topright",
           legend = c(legend_metric, "GEMMA all", "GEMMA male", "Campbell genes"),
           col = c(legend_metric_col, COL_GEMMA_ALL, COL_GEMMA_MALE, COL_CAMPBELL),
           pch = c(rep(20, length(legend_metric)), 24, 23, NA),
           pt.bg = c(rep(NA, length(legend_metric)), COL_GEMMA_ALL, COL_GEMMA_MALE, NA),
           lty = c(rep(NA, length(legend_metric) + 2), 3),
           bty = "n", cex = 0.75, pt.cex = 0.9)
  }

  png(png_file, width = MAIN_WIDTH, height = MAIN_HEIGHT, res = PNG_RES)
  plot_fun()
  dev.off()

  pdf(pdf_file, width = 20, height = 7.4)
  plot_fun()
  dev.off()

  thresholds_out <- data.table(
    figure = outfile_prefix,
    threshold = c("q001", "q01", "q99", "q999"),
    value = c(q001, q01, q99, q999)
  )

  list(png = png_file, pdf = pdf_file, thresholds = thresholds_out)
}

# ============================================================
# Read and standardize scan tables
# ============================================================

msg("Reading scan tables")
fst_raw <- fread(FST_50KB)
ho_raw  <- fread(HO_50KB)
pi_raw  <- fread(PI_50KB)
taj_raw <- fread(TAJ_50KB)

# FST columns
fst_chrom <- first_existing_col(fst_raw, c("CHR", "CHROM"))
fst_start <- first_existing_col(fst_raw, c("window_start", "BIN_START", "START"))
fst_end   <- first_existing_col(fst_raw, c("window_end", "BIN_END", "END"))
fst_value <- first_existing_col(fst_raw, c("Mean_FST", "mean_FST", "FST"))
fst_snps  <- first_existing_col(fst_raw, c("N_SNPs", "N_VARIANTS"))
fst <- standardize_windows(fst_raw, fst_chrom, fst_start, fst_end, fst_value, fst_snps, min_snps = 2)

# Ho columns
ho_chrom <- first_existing_col(ho_raw, c("CHROM", "CHR"))
ho_start <- first_existing_col(ho_raw, c("BIN_START", "window_start", "START"))
ho_end   <- first_existing_col(ho_raw, c("BIN_END", "window_end", "END"))
ho_value <- first_existing_col(ho_raw, c("delta_Ho_TX_minus_TT", "delta_Ho", "Ho_TX_minus_TT"))
ho_snps  <- first_existing_col(ho_raw, c("N_SNPs", "N_VARIANTS"))
ho <- standardize_windows(ho_raw, ho_chrom, ho_start, ho_end, ho_value, ho_snps, min_snps = 10)

# pi columns
pi_chrom <- first_existing_col(pi_raw, c("CHROM", "CHR"))
pi_start <- first_existing_col(pi_raw, c("BIN_START", "window_start", "START"))
pi_end   <- first_existing_col(pi_raw, c("BIN_END", "window_end", "END"))
pi_value <- first_existing_col(pi_raw, c("delta_pi", "delta_PI_TX_minus_TT", "PI_TX_minus_TT"))
pi <- standardize_windows(pi_raw, pi_chrom, pi_start, pi_end, pi_value)

# TajimaD columns
taj_chrom <- first_existing_col(taj_raw, c("CHROM", "CHR"))
taj_start <- first_existing_col(taj_raw, c("BIN_START", "window_start", "START"))
taj_end   <- first_existing_col(taj_raw, c("BIN_END", "window_end", "END"))
taj_value <- find_tajima_value_col(taj_raw)
taj <- standardize_windows(taj_raw, taj_chrom, taj_start, taj_end, taj_value)

msg("Windows after filtering:",
    "FST", nrow(fst),
    "Ho", nrow(ho),
    "pi", nrow(pi),
    "TajimaD", nrow(taj))

layout <- make_layout(fst, ho, pi, taj)
fst <- add_xpos_windows(fst, layout)
ho  <- add_xpos_windows(ho, layout)
pi  <- add_xpos_windows(pi, layout)
taj <- add_xpos_windows(taj, layout)

# ============================================================
# Markers: GEMMA and Campbell genes
# ============================================================

gemma_all <- read_gemma_hits(GEMMA_ALL, "GEMMA_all_PC3")
gemma_male <- read_gemma_hits(GEMMA_MALE, "GEMMA_male_PC3")

gemma_all[, marker_type := "GEMMA_all_Bonferroni"]
gemma_male[, marker_type := "GEMMA_male_Bonferroni"]

gemma_all <- attach_layout_to_markers(gemma_all, layout)
gemma_male <- attach_layout_to_markers(gemma_male, layout)

campbell <- read_campbell_from_gff(GFF, MAP)
fallback <- fallback_campbell()
if (nrow(campbell) == 0) {
  msg("No Campbell genes detected from GFF; using fallback positions for available genes")
  campbell <- fallback
} else {
  missing_symbols <- setdiff(fallback$Campbell_gene, campbell$Campbell_gene)
  if (length(missing_symbols) > 0) {
    campbell <- rbindlist(list(campbell, fallback[Campbell_gene %in% missing_symbols]), fill = TRUE)
  }
}
campbell <- attach_layout_to_markers(campbell, layout)

# Marker table
marker_table <- rbindlist(
  list(
    gemma_all[, .(marker_type, model, SNP_ID, CHROM, POS, P, chr_num, xpos, short_label,
                  Campbell_gene = NA_character_, gene_start = NA_real_, gene_end = NA_real_)],
    gemma_male[, .(marker_type, model, SNP_ID, CHROM, POS, P, chr_num, xpos, short_label,
                   Campbell_gene = NA_character_, gene_start = NA_real_, gene_end = NA_real_)],
    campbell[, .(marker_type, model, SNP_ID, CHROM, POS, P, chr_num, xpos, short_label,
                 Campbell_gene, gene_start = START, gene_end = END)]
  ),
  fill = TRUE
)

fwrite(marker_table,
       file.path(TABDIR, "Genomewide_scan_marker_positions_GEMMA_Campbell.tsv"),
       sep = "\t")

msg("Markers:",
    "GEMMA_all", nrow(gemma_all),
    "GEMMA_male", nrow(gemma_male),
    "Campbell", nrow(campbell))

# ============================================================
# Plot all scans
# ============================================================

res <- list()
res[["FST"]] <- plot_scan(
  fst, layout, gemma_all, gemma_male, campbell,
  outfile_prefix = "Fig1_FST_50kb_GEMMA_Campbell_annotated",
  main = "Genome-wide FST scan TT vs TX, 50 kb windows",
  ylab = expression("Mean " * F[ST]),
  mode = "high"
)

res[["Ho"]] <- plot_scan(
  ho, layout, gemma_all, gemma_male, campbell,
  outfile_prefix = "Fig2_delta_Ho_50kb_GEMMA_Campbell_annotated",
  main = "Genome-wide observed heterozygosity difference, 50 kb windows",
  ylab = expression(Delta~Ho~"(TX - TT)"),
  mode = "low"
)

res[["pi"]] <- plot_scan(
  pi, layout, gemma_all, gemma_male, campbell,
  outfile_prefix = "Fig3_delta_pi_50kb_GEMMA_Campbell_annotated",
  main = "Genome-wide nucleotide diversity difference, 50 kb windows",
  ylab = expression(Delta*pi~"(TX - TT)"),
  mode = "two_sided"
)

res[["TajimaD"]] <- plot_scan(
  taj, layout, gemma_all, gemma_male, campbell,
  outfile_prefix = "Fig4_delta_TajimaD_50kb_GEMMA_Campbell_annotated",
  main = "Genome-wide Tajima's D difference, 50 kb windows",
  ylab = expression(Delta~"Tajima's D"~"(TX - TT)"),
  mode = "two_sided"
)

thresholds <- rbindlist(lapply(names(res), function(nm) {
  x <- res[[nm]]$thresholds
  x[, scan := nm]
  x
}), fill = TRUE)

fwrite(thresholds,
       file.path(TABDIR, "Genomewide_scan_thresholds_GEMMA_Campbell_annotated.tsv"),
       sep = "\t")

summary <- data.table(
  item = c(
    "FST_windows", "Ho_windows", "pi_windows", "TajimaD_windows",
    "GEMMA_all_markers", "GEMMA_male_markers", "Campbell_gene_markers",
    "Campbell_genes", "Output_directory"
  ),
  value = c(
    nrow(fst), nrow(ho), nrow(pi), nrow(taj),
    nrow(gemma_all), nrow(gemma_male), nrow(campbell),
    paste(campbell$Campbell_gene, collapse = ";"),
    OUTDIR
  )
)

fwrite(summary,
       file.path(TABDIR, "Genomewide_scan_GEMMA_Campbell_annotated_summary.tsv"),
       sep = "\t")

cat("\nDone. Annotated figures written to:\n", OUTDIR, "\n", sep = "")
cat("Marker table:\n", file.path(TABDIR, "Genomewide_scan_marker_positions_GEMMA_Campbell.tsv"), "\n", sep = "")
cat("Threshold table:\n", file.path(TABDIR, "Genomewide_scan_thresholds_GEMMA_Campbell_annotated.tsv"), "\n", sep = "")
