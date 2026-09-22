suppressPackageStartupMessages(library(data.table))

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
SNPDIR <- file.path(CAMP, "snp_level_campbell")
OUTDIR <- file.path(CAMP, "strict_HetDev_venn")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

fread_any <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)))
  } else {
    fread(path)
  }
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

make_pattern <- function(f, d, m, h) {
  labs <- c("FST", "DXY", "muLD", "HetDev_TX")
  vals <- c(f, d, m, h)
  x <- labs[vals]
  if (length(x) == 0) "none" else paste(x, collapse = "&")
}

summarize_one <- function(file, label) {
  cat("\nReading:", file, "\n")
  mem <- fread_any(file)

  mem[, FST_b := as_bool(FST)]
  mem[, DXY_b := as_bool(DXY)]
  mem[, muLD_b := as_bool(muLD)]
  mem[, HetDev_projected_b := as_bool(HetDev_TX)]
  mem[, Campbell_Het_candidate_b := as_bool(Campbell_Het_candidate)]

  # Version stricte :
  # un SNP est HetDev seulement s'il est dans une fenêtre HetDev
  # ET s'il satisfait lui-même le critère Campbell_Het_candidate.
  mem[, HetDev_strict_b := HetDev_projected_b & Campbell_Het_candidate_b]

  for (mode in c("projected", "strictHet")) {
    h <- if (mode == "projected") {
      mem$HetDev_projected_b
    } else {
      mem$HetDev_strict_b
    }

    tmp <- data.table(
      snp_id = mem$snp_id,
      CHROM = mem$CHROM,
      POS = as.numeric(mem$POS),
      FST = mem$FST_b,
      DXY = mem$DXY_b,
      muLD = mem$muLD_b,
      HetDev_TX = h,
      HetDev_projected = mem$HetDev_projected_b,
      Campbell_Het_candidate = mem$Campbell_Het_candidate_b
    )

    tmp[, n_metrics := FST + DXY + muLD + HetDev_TX]
    tmp[, pattern := mapply(make_pattern, FST, DXY, muLD, HetDev_TX)]

    tmp_sel <- tmp[n_metrics > 0]
    counts <- tmp_sel[, .N, by = .(pattern, n_metrics)][order(-n_metrics, pattern)]

    fwrite(
      tmp_sel,
      file.path(OUTDIR, paste0("membership_", label, "_", mode, ".tsv")),
      sep = "\t"
    )

    fwrite(
      counts,
      file.path(OUTDIR, paste0("pattern_counts_", label, "_", mode, ".tsv")),
      sep = "\t"
    )

    fwrite(
      tmp_sel[n_metrics >= 3 & HetDev_TX == TRUE],
      file.path(OUTDIR, paste0("SNPs_at_least3_including_HetDev_", label, "_", mode, ".tsv")),
      sep = "\t"
    )

    cat("\n===", label, mode, "===\n")
    print(counts)
    cat("SNPs >=3 metrics including HetDev:", nrow(tmp_sel[n_metrics >= 3 & HetDev_TX == TRUE]), "\n")
  }
}

summarize_one(
  file.path(SNPDIR, "Campbell_SNP_membership_1pct_Campbell_threshold.tsv.gz"),
  "1pct"
)

summarize_one(
  file.path(SNPDIR, "Campbell_SNP_membership_5pct_Campbell.tsv.gz"),
  "5pct"
)
