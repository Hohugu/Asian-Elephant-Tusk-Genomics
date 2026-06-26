#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST"

infile <- file.path(
  base,
  "tables",
  "TT_vs_TX_FST_50kb_genomewide_highlight_candidate_windows_p999.tsv"
)

outfile <- file.path(
  base,
  "tables",
  "TT_vs_TX_FST_50kb_candidate_regions_merged.tsv"
)

gap_bp <- 100000

x <- fread(infile)

setorder(x, CHR, window_start)

x[, region_id := NA_integer_]

current_region <- 1
x$region_id[1] <- current_region

for (i in 2:nrow(x)) {

  same_chr <- x$CHR[i] == x$CHR[i - 1]

  gap <- x$window_start[i] - x$window_end[i - 1]

  if (same_chr && gap <= gap_bp) {
    x$region_id[i] <- current_region
  } else {
    current_region <- current_region + 1
    x$region_id[i] <- current_region
  }
}

regions <- x[
  ,
  .(
    CHR = CHR[1],
    Start = min(window_start),
    End = max(window_end),
    Length_bp = max(window_end) - min(window_start) + 1,
    N_windows = .N,
    N_SNPs = sum(N_SNPs),
    Mean_FST = mean(Mean_FST, na.rm = TRUE),
    Max_Mean_FST = max(Mean_FST, na.rm = TRUE),
    Max_SNP_FST = max(Max_FST, na.rm = TRUE),
    Top_SNP = Top_SNP[which.max(Max_FST)],
    Top_SNP_POS = Top_SNP_POS[which.max(Max_FST)]
  ),
  by = region_id
]

setorder(regions, -Max_Mean_FST)

regions[, Region_ID := paste0("FST_region_", .I)]

setcolorder(
  regions,
  c(
    "Region_ID",
    "CHR",
    "Start",
    "End",
    "Length_bp",
    "N_windows",
    "N_SNPs",
    "Mean_FST",
    "Max_Mean_FST",
    "Max_SNP_FST",
    "Top_SNP",
    "Top_SNP_POS"
  )
)

regions[, Start := as.integer(Start)]
regions[, End := as.integer(End)]
regions[, Length_bp := as.integer(Length_bp)]
regions[, Top_SNP_POS := as.integer(Top_SNP_POS)]
regions[, region_id := NULL]

fwrite(regions, outfile, sep = "\t")

cat("Done\n")
cat("Input p999 windows:", nrow(x), "\n")
cat("Merged candidate regions:", nrow(regions), "\n")
cat("Output:", outfile, "\n")
