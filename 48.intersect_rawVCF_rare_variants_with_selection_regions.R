#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

rare_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_differentiated_variants.tsv"
)

conv_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv"
)

xy_file <- file.path(
  base,
  "Genetics_Analysis/X_chromosome/XY_selection_signals_50kb_integrated.tsv"
)

gwas_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/ALL_GWAS_Bonferroni_SNPs_annotated_fullGFF.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/rawVCF_rare_variants_overlap_selection_GWAS.tsv"
)

summary_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/rawVCF_rare_variants_overlap_selection_GWAS_summary.tsv"
)

rare <- fread(rare_file)
conv <- fread(conv_file)
xy <- fread(xy_file)
gwas <- fread(gwas_file)

conv_regions <- conv[, .(
  region_type = "convergent_region",
  region_id = Region_ID,
  CHROM,
  START,
  END
)]

xy_regions <- xy[Evidence_score >= 3, .(
  region_type = "XY_score3_region",
  region_id = paste0(CHROM, ":", BIN_START, "-", BIN_END),
  CHROM,
  START = BIN_START,
  END = BIN_END
)]

regions <- rbindlist(list(conv_regions, xy_regions), fill = TRUE)

rare[, in_convergent_region := FALSE]
rare[, convergent_region_id := NA_character_]
rare[, in_XY_score3_region := FALSE]
rare[, XY_score3_region_id := NA_character_]

for (i in seq_len(nrow(regions))) {
  r <- regions[i]

  idx <- rare$CHROM == r$CHROM &
    rare$POS >= r$START &
    rare$POS <= r$END

  if (any(idx)) {
    if (r$region_type == "convergent_region") {
      rare[idx, in_convergent_region := TRUE]
      rare[idx, convergent_region_id := r$region_id]
    }

    if (r$region_type == "XY_score3_region") {
      rare[idx, in_XY_score3_region := TRUE]
      rare[idx, XY_score3_region_id := r$region_id]
    }
  }
}

gwas_snps <- unique(gwas$SNP)

rare[, is_GWAS_Bonferroni_SNP := SNP %in% gwas_snps]

summary <- data.table(
  metric = c(
    "rare_differentiated_total",
    "in_convergent_region",
    "in_XY_score3_region",
    "is_GWAS_Bonferroni_SNP",
    "in_convergent_and_GWAS",
    "in_XY_and_GWAS"
  ),
  N = c(
    nrow(rare),
    sum(rare$in_convergent_region),
    sum(rare$in_XY_score3_region),
    sum(rare$is_GWAS_Bonferroni_SNP),
    sum(rare$in_convergent_region & rare$is_GWAS_Bonferroni_SNP),
    sum(rare$in_XY_score3_region & rare$is_GWAS_Bonferroni_SNP)
  )
)

setorder(rare, -abs_delta_AF)

fwrite(rare, out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(summary)
