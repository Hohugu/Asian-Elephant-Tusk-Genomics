#!/usr/bin/env Rscript

camp <- Sys.getenv("CAMP")
if (camp == "") stop("CAMP not set")

outdir <- file.path(camp, "snp_level_campbell")
lib <- file.path(camp, "Rlibs")
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", lib = lib, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", lib = lib, repos = "https://cloud.r-project.org")
}

library(data.table)
library(VennDiagram)
library(grid)

metrics <- c("FST", "DXY", "muLD", "HetDev_TX")

get_category_count <- function(dt, category) {
  x <- dt[category == !!category, count]
  if (length(x) == 0) return(0)
  as.numeric(x[1])
}

contains_all <- function(pattern, combo) {
  parts <- unlist(strsplit(pattern, "&", fixed = TRUE))
  all(combo %in% parts)
}

make_venn <- function(label, title_suffix) {
  counts_file <- file.path(outdir, paste0("Campbell_SNP_counts_", label, ".tsv"))
  pair_file <- file.path(outdir, paste0("Campbell_SNP_pairwise_counts_", label, ".tsv"))
  pattern_file <- file.path(outdir, paste0("Campbell_SNP_pattern_counts_", label, ".tsv"))

  counts <- fread(counts_file)
  pairs <- fread(pair_file)
  patterns <- fread(pattern_file)

  get_count <- function(cat) {
    x <- counts[category == cat, count]
    if (length(x) == 0) return(0)
    as.numeric(x[1])
  }

  get_pair <- function(a, b) {
    x <- pairs[metric1 == a & metric2 == b, count]
    if (length(x) == 0) x <- pairs[metric1 == b & metric2 == a, count]
    if (length(x) == 0) return(0)
    as.numeric(x[1])
  }

  get_combo <- function(combo) {
    sum(patterns[sapply(pattern, contains_all, combo = combo), count])
  }

  area1 <- get_count("FST_SNPs")
  area2 <- get_count("DXY_SNPs")
  area3 <- get_count("muLD_SNPs")
  area4 <- get_count("HetDev_TX_SNPs")

  n12 <- get_pair("FST", "DXY")
  n13 <- get_pair("FST", "muLD")
  n14 <- get_pair("FST", "HetDev_TX")
  n23 <- get_pair("DXY", "muLD")
  n24 <- get_pair("DXY", "HetDev_TX")
  n34 <- get_pair("muLD", "HetDev_TX")

  n123 <- get_combo(c("FST", "DXY", "muLD"))
  n124 <- get_combo(c("FST", "DXY", "HetDev_TX"))
  n134 <- get_combo(c("FST", "muLD", "HetDev_TX"))
  n234 <- get_combo(c("DXY", "muLD", "HetDev_TX"))
  n1234 <- get_combo(c("FST", "DXY", "muLD", "HetDev_TX"))

  summary <- data.table(
    threshold = label,
    category = c(
      "FST", "DXY", "muLD", "HetDev_TX",
      "FST&DXY", "FST&muLD", "FST&HetDev_TX",
      "DXY&muLD", "DXY&HetDev_TX", "muLD&HetDev_TX",
      "FST&DXY&muLD", "FST&DXY&HetDev_TX",
      "FST&muLD&HetDev_TX", "DXY&muLD&HetDev_TX",
      "FST&DXY&muLD&HetDev_TX"
    ),
    count = c(
      area1, area2, area3, area4,
      n12, n13, n14,
      n23, n24, n34,
      n123, n124, n134, n234,
      n1234
    )
  )

  fwrite(
    summary,
    file.path(outdir, paste0("Campbell_SNP_level_venn_counts_", label, ".tsv")),
    sep = "\t"
  )

  png_file <- file.path(outdir, paste0("Campbell_SNP_level_Venn_", label, ".png"))

  png(png_file, width = 3200, height = 3200, res = 300)
  grid.newpage()

  venn <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = c("FST", "DXY", "muLD", "HetDev_TX"),
    fill = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
    alpha = 0.35,
    cex = 0.75,
    cat.cex = 0.9,
    margin = 0.08,
    main = paste0("Campbell-style SNP-level Venn, ", title_suffix)
  )

  grid.draw(venn)
  dev.off()

  message("Wrote: ", png_file)
}

make_venn("5pct_Campbell", "top 5% tails")
make_venn("1pct_Campbell_threshold", "top 1% tails")
