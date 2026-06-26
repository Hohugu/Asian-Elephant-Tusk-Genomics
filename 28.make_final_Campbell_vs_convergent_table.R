#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

conv_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/final_candidate_regions_summary.tsv"
)

camp_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/Campbell_candidates_selection_signals.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/final_Campbell_vs_convergent_selection_summary.tsv"
)

conv <- fread(conv_file)
camp <- fread(camp_file)

conv2 <- conv[
  ,
  .(
    Category = "Convergent_region",
    Name = Region_ID,
    Gene_or_region = fifelse(
      !is.na(Protein_coding_descriptions) & Protein_coding_descriptions != "",
      Protein_coding_descriptions,
      Gene_description
    ),
    CHROM,
    START,
    END,
    Max_FST = Max_SNP_FST,
    delta_pi,
    delta_TajimaD,
    delta_Ho = delta_Ho_TX_minus_TT,
    mean_r2_TT,
    mean_r2_TX,
    delta_LD = delta_mean_r2_TX_minus_TT,
    Top_SNP = Top_deltaAF_SNP,
    AF_TT,
    AF_TX,
    delta_AF = delta_AF_TX_minus_TT,
    Evidence_score = Evidence_score_with_Het
  )
]

camp2 <- camp[
  ,
  .(
    Category = "Campbell_candidate",
    Name = Campbell_name,
    Gene_or_region = description,
    CHROM = CM_chr,
    START = start,
    END = end,
    Max_FST,
    delta_pi,
    delta_TajimaD,
    delta_Ho = delta_Ho_TX_minus_TT,
    mean_r2_TT,
    mean_r2_TX,
    delta_LD = delta_mean_r2_TX_minus_TT,
    Top_SNP = NA_character_,
    AF_TT = NA_real_,
    AF_TX = NA_real_,
    delta_AF = NA_real_,
    Evidence_score = NA_integer_
  )
]

res <- rbindlist(list(conv2, camp2), fill = TRUE)

res[, Signal_summary := paste0(
  fifelse(!is.na(Max_FST) & Max_FST >= 0.10, "FST_high;", ""),
  fifelse(!is.na(delta_pi) & delta_pi < 0, "pi_TX_lower;", ""),
  fifelse(!is.na(delta_TajimaD) & delta_TajimaD < 0, "TajimaD_TX_lower;", ""),
  fifelse(!is.na(delta_Ho) & delta_Ho < 0, "Ho_TX_lower;", ""),
  fifelse(!is.na(delta_LD) & delta_LD > 0, "LD_TX_higher;", "")
)]

res[, Campbell_like_score := 0L]
res[!is.na(Max_FST) & Max_FST >= 0.10, Campbell_like_score := Campbell_like_score + 1L]
res[!is.na(delta_pi) & delta_pi < 0, Campbell_like_score := Campbell_like_score + 1L]
res[!is.na(delta_TajimaD) & delta_TajimaD < 0, Campbell_like_score := Campbell_like_score + 1L]
res[!is.na(delta_Ho) & delta_Ho < 0, Campbell_like_score := Campbell_like_score + 1L]
res[!is.na(delta_LD) & delta_LD > 0, Campbell_like_score := Campbell_like_score + 1L]

setorder(res, -Campbell_like_score, Category, Name)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(res)
