#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

candidate_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD.tsv"
)

het_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv"
)

cand <- fread(candidate_file)
het <- fread(het_file)

cand[, BIN_START := START]
cand[, BIN_END := END]

het_sub <- het[
  ,
  .(
    CHROM,
    BIN_START,
    BIN_END,
    Ho_TT,
    Ho_TX,
    delta_Ho_TX_minus_TT
  )
]

res <- merge(
  cand,
  het_sub,
  by = c("CHROM", "BIN_START", "BIN_END"),
  all.x = TRUE
)

res[, Het_direction := fifelse(
  delta_Ho_TX_minus_TT < 0,
  "TX lower heterozygosity",
  "TX higher heterozygosity"
)]

res[, Evidence_score_with_Het := Evidence_score]
res[delta_Ho_TX_minus_TT < 0, Evidence_score_with_Het := Evidence_score_with_Het + 1L]

setorder(res, -Evidence_score_with_Het, delta_pi)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(res[, .(
  Region_ID, CHROM, START, END,
  Max_SNP_FST,
  delta_pi,
  delta_TajimaD,
  Ho_TT,
  Ho_TX,
  delta_Ho_TX_minus_TT,
  Het_direction,
  Evidence_score_with_Het,
  Protein_coding_descriptions
)])
