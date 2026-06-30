#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

outdir <- file.path(base, "Genetics_Analysis/Functional_annotation/tables")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# GWAS Bonferroni SNPs
# -------------------------

gwas_all_file <- file.path(
  base,
  "GWAS/result/GEMMA_PC3/output/GEMMA_PC3_Bonferroni_hits.tsv"
)

gwas_male_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/GEMMA_male_Bonferroni_hits.tsv"
)

gwas_all <- fread(gwas_all_file)
gwas_male <- fread(gwas_male_file)

gwas_all <- gwas_all[, .(
  GWAS_model = "GEMMA_all_PC3",
  CHROM = chr,
  POS = ps,
  SNP = rs,
  P = p_wald
)]

gwas_male <- gwas_male[, .(
  GWAS_model = "GEMMA_male_PC3",
  CHROM = chr,
  POS = ps,
  SNP = rs,
  P = p_wald
)]

gwas <- rbindlist(list(gwas_all, gwas_male), fill = TRUE)

# -------------------------
# Selection regions/windows
# -------------------------

fst_regions_file <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_candidate_regions_merged.tsv"
)

convergent_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv"
)

xy_file <- file.path(
  base,
  "Genetics_Analysis/X_chromosome/XY_selection_signals_50kb_integrated.tsv"
)

pi_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_pi_50kb_top100_reduced_TX.tsv"
)

taj_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_top100_more_negative_TX.tsv"
)

ho_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_delta_Ho_50kb_top100_reduced_TX.tsv"
)

fst_regions <- fread(fst_regions_file)
conv <- fread(convergent_file)
xy <- fread(xy_file)
pi_top <- fread(pi_file)
taj_top <- fread(taj_file)
ho_top <- fread(ho_file)

# Standardize coordinates
fst_regions_std <- fst_regions[, .(
  Region_type = "FST_candidate_region",
  Region_ID = Region_ID,
  CHROM = CHR,
  START = as.integer(Start),
  END = as.integer(End)
)]

conv_std <- conv[, .(
  Region_type = "Convergent_FST_pi_TajimaD_Ho",
  Region_ID,
  CHROM,
  START,
  END
)]

xy_std <- xy[Evidence_score >= 3, .(
  Region_type = "XY_selection_score3",
  Region_ID = paste0(CHROM, ":", BIN_START, "-", BIN_END),
  CHROM,
  START = BIN_START,
  END = BIN_END
)]

pi_std <- pi_top[, .(
  Region_type = "pi_top100_reduced_TX",
  Region_ID = paste0(CHROM, ":", BIN_START, "-", BIN_END),
  CHROM,
  START = BIN_START,
  END = BIN_END
)]

taj_std <- taj_top[, .(
  Region_type = "TajimaD_top100_more_negative_TX",
  Region_ID = paste0(CHROM, ":", BIN_START, "-", BIN_END),
  CHROM,
  START = BIN_START,
  END = BIN_END
)]

ho_std <- ho_top[, .(
  Region_type = "Ho_top100_reduced_TX",
  Region_ID = paste0(CHROM, ":", BIN_START, "-", BIN_END),
  CHROM,
  START = BIN_START,
  END = BIN_END
)]

regions <- rbindlist(
  list(fst_regions_std, conv_std, xy_std, pi_std, taj_std, ho_std),
  fill = TRUE
)

regions <- unique(regions)

# -------------------------
# Overlap SNPs with regions
# -------------------------

overlap_one <- function(snp_dt, region_dt, type_name) {
  r <- region_dt[Region_type == type_name]

  if (nrow(r) == 0) {
    return(snp_dt[, .(
      GWAS_model, SNP,
      hit = FALSE,
      region_ids = NA_character_
    )])
  }

  out <- snp_dt[
    ,
    {
      hits <- r[CHROM == .BY$CHROM & START <= .BY$POS & END >= .BY$POS]
      .(
        hit = nrow(hits) > 0,
        region_ids = ifelse(nrow(hits) > 0, paste(unique(hits$Region_ID), collapse = ";"), NA_character_)
      )
    },
    by = .(GWAS_model, SNP, CHROM, POS)
  ]

  setnames(out, c("hit", "region_ids"), c(paste0("in_", type_name), paste0(type_name, "_IDs")))
  out[, .SD, .SDcols = c("GWAS_model", "SNP", paste0("in_", type_name), paste0(type_name, "_IDs"))]
}

types <- unique(regions$Region_type)

annot <- copy(gwas)

for (tp in types) {
  ov <- overlap_one(gwas, regions, tp)
  annot <- merge(annot, ov, by = c("GWAS_model", "SNP"), all.x = TRUE)
}

# Any selection overlap
hit_cols <- grep("^in_", names(annot), value = TRUE)
annot[, Any_selection_overlap := rowSums(.SD == TRUE, na.rm = TRUE) > 0, .SDcols = hit_cols]
annot[, N_selection_categories := rowSums(.SD == TRUE, na.rm = TRUE), .SDcols = hit_cols]

# -------------------------
# Nearest selection region distance
# -------------------------

nearest_list <- list()

for (i in seq_len(nrow(gwas))) {
  s <- gwas[i]
  same_chr <- regions[CHROM == s$CHROM]

  if (nrow(same_chr) == 0) {
    nearest_list[[i]] <- data.table(
      GWAS_model = s$GWAS_model,
      SNP = s$SNP,
      nearest_selection_region = NA_character_,
      nearest_selection_type = NA_character_,
      distance_to_nearest_selection_bp = NA_real_
    )
  } else {
    same_chr[, dist := pmax(START - s$POS, s$POS - END, 0)]
    nearest <- same_chr[which.min(dist)]
    nearest_list[[i]] <- data.table(
      GWAS_model = s$GWAS_model,
      SNP = s$SNP,
      nearest_selection_region = nearest$Region_ID,
      nearest_selection_type = nearest$Region_type,
      distance_to_nearest_selection_bp = nearest$dist
    )
  }
}

nearest <- rbindlist(nearest_list)

annot <- merge(annot, nearest, by = c("GWAS_model", "SNP"), all.x = TRUE)

setorder(annot, GWAS_model, P)

# -------------------------
# Summary
# -------------------------

summary <- annot[
  ,
  c(
    .(
      N_Bonferroni_SNPs = .N,
      N_any_selection_overlap = sum(Any_selection_overlap),
      Percent_any_selection_overlap = round(100 * mean(Any_selection_overlap), 2)
    ),
    lapply(.SD, function(x) sum(x == TRUE, na.rm = TRUE))
  ),
  by = GWAS_model,
  .SDcols = hit_cols
]

# -------------------------
# Outputs
# -------------------------

out_snp <- file.path(outdir, "GEMMA_Bonferroni_all_and_male_vs_selection_SNP_level.tsv")
out_sum <- file.path(outdir, "GEMMA_Bonferroni_all_and_male_vs_selection_summary.tsv")

fwrite(annot, out_snp, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, out_sum, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("SNP-level:", out_snp, "\n")
cat("Summary:", out_sum, "\n\n")

print(summary)
