#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

effects_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_CDS_variant_effects.tsv"
)

rare_annot_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_regulatory_priority_annotation.tsv"
)

out_unique <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_CDS_missense_unique_variants.tsv"
)

out_high <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_CDS_missense_high_priority.tsv"
)

out_summary <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_CDS_missense_unique_summary.tsv"
)

x <- fread(effects_file)
rare <- fread(rare_annot_file)

miss <- x[Predicted_effect == "missense"]

miss_unique <- miss[, .(
  rare_class = first(rare_class),
  AF_TT = first(AF_TT),
  AF_TX = first(AF_TX),
  abs_delta_AF = first(abs_delta_AF),
  Genes = paste(sort(unique(Gene_ID)), collapse = ";"),
  Products = paste(sort(unique(Product)), collapse = ";"),
  AA_changes = paste(sort(unique(paste0(AA_REF, AA_position, AA_ALT))), collapse = ";"),
  Transcripts = paste(sort(unique(Transcript)), collapse = ";"),
  Protein_IDs = paste(sort(unique(Protein_ID)), collapse = ";"),
  N_transcript_rows = .N,
  REF_match_values = paste(sort(unique(REF_match)), collapse = ";"),
  nearest_gene_TSS_from_effects = paste(sort(unique(nearest_gene_TSS)), collapse = ";"),
  nearest_gene_description_from_effects = paste(sort(unique(nearest_gene_description_TSS)), collapse = ";"),
  rare_variant_priority_reason_from_effects = paste(sort(unique(rare_variant_priority_reason)), collapse = ";")
), by = .(SNP, CHROM, POS, REF, ALT)]

rare_keep_cols <- intersect(
  c(
    "SNP",
    "regulatory_class",
    "in_convergent_region",
    "convergent_region_id",
    "in_XY_score3_region",
    "XY_score3_region_id",
    "nearest_gene_TSS",
    "nearest_gene_description_TSS",
    "nearest_gene_biotype_TSS",
    "distance_to_nearest_TSS_bp",
    "distance_to_gene_body_bp",
    "SNP_location_relative_to_gene",
    "nearest_gene_is_GWAS_priority",
    "rare_variant_priority_reason"
  ),
  names(rare)
)

rare_small <- rare[, ..rare_keep_cols]

miss_unique <- merge(
  miss_unique,
  rare_small,
  by = "SNP",
  all.x = TRUE
)

bio_pattern <- paste(
  c(
    "tooth",
    "teeth",
    "odont",
    "dentin",
    "enamel",
    "ameloblast",
    "cranio",
    "facial",
    "neural crest",
    "bone",
    "osteo",
    "cartilage",
    "BMP",
    "WNT",
    "SHH",
    "FGF",
    "TGF",
    "RUNX",
    "DLX",
    "MSX",
    "PAX",
    "EDA",
    "EDAR",
    "AXIN",
    "LRP",
    "TCF",
    "collagen",
    "matrix",
    "morphogen",
    "development",
    "growth",
    "polycystic",
    "fibrocystin",
    "rho GTPase",
    "GTPase",
    "teashirt",
    "ubiquitin",
    "histone",
    "methyltransferase",
    "carbonyl reductase",
    "EGF"
  ),
  collapse = "|"
)

miss_unique[, text_for_keyword := paste(
  Genes,
  Products,
  nearest_gene_TSS,
  nearest_gene_description_TSS,
  sep = " "
)]

miss_unique[, has_biological_keyword := grepl(
  bio_pattern,
  text_for_keyword,
  ignore.case = TRUE
)]

miss_unique[, missense_priority_score := 0L]

miss_unique[REF_match_values == "TRUE",
  missense_priority_score := missense_priority_score + 1L]

miss_unique[distance_to_nearest_TSS_bp <= 2000,
  missense_priority_score := missense_priority_score + 2L]

miss_unique[distance_to_nearest_TSS_bp > 2000 &
              distance_to_nearest_TSS_bp <= 10000,
  missense_priority_score := missense_priority_score + 1L]

miss_unique[nearest_gene_is_GWAS_priority == TRUE,
  missense_priority_score := missense_priority_score + 3L]

miss_unique[in_convergent_region == TRUE,
  missense_priority_score := missense_priority_score + 3L]

miss_unique[in_XY_score3_region == TRUE,
  missense_priority_score := missense_priority_score + 2L]

miss_unique[has_biological_keyword == TRUE,
  missense_priority_score := missense_priority_score + 2L]

miss_unique[, missense_priority_reason := ""]

miss_unique[REF_match_values == "TRUE",
  missense_priority_reason := paste0(missense_priority_reason, "REF_match_TRUE;")]

miss_unique[distance_to_nearest_TSS_bp <= 2000,
  missense_priority_reason := paste0(missense_priority_reason, "TSS_within_2kb;")]

miss_unique[distance_to_nearest_TSS_bp > 2000 &
              distance_to_nearest_TSS_bp <= 10000,
  missense_priority_reason := paste0(missense_priority_reason, "TSS_within_10kb;")]

miss_unique[nearest_gene_is_GWAS_priority == TRUE,
  missense_priority_reason := paste0(missense_priority_reason, "GWAS_priority_gene;")]

miss_unique[in_convergent_region == TRUE,
  missense_priority_reason := paste0(missense_priority_reason, "convergent_selection;")]

miss_unique[in_XY_score3_region == TRUE,
  missense_priority_reason := paste0(missense_priority_reason, "XY_selection;")]

miss_unique[has_biological_keyword == TRUE,
  missense_priority_reason := paste0(missense_priority_reason, "biological_keyword;")]

miss_unique[, missense_priority_reason := sub(";$", "", missense_priority_reason)]

miss_unique[missense_priority_reason == "",
  missense_priority_reason := "missense_only"]

setorder(
  miss_unique,
  -missense_priority_score,
  -abs_delta_AF,
  CHROM,
  POS
)

high <- miss_unique[
  missense_priority_score >= 4 |
    nearest_gene_is_GWAS_priority == TRUE |
    in_convergent_region == TRUE |
    in_XY_score3_region == TRUE
]

summary <- data.table(
  metric = c(
    "unique_missense_variants",
    "high_priority_missense_variants",
    "REF_match_TRUE",
    "within_2kb_TSS",
    "within_10kb_TSS",
    "nearest_GWAS_priority_gene",
    "in_convergent_region",
    "in_XY_score3_region",
    "has_biological_keyword"
  ),
  N = c(
    nrow(miss_unique),
    nrow(high),
    sum(miss_unique$REF_match_values == "TRUE", na.rm = TRUE),
    sum(miss_unique$distance_to_nearest_TSS_bp <= 2000, na.rm = TRUE),
    sum(miss_unique$distance_to_nearest_TSS_bp <= 10000, na.rm = TRUE),
    sum(miss_unique$nearest_gene_is_GWAS_priority == TRUE, na.rm = TRUE),
    sum(miss_unique$in_convergent_region == TRUE, na.rm = TRUE),
    sum(miss_unique$in_XY_score3_region == TRUE, na.rm = TRUE),
    sum(miss_unique$has_biological_keyword == TRUE, na.rm = TRUE)
  )
)

fwrite(miss_unique, out_unique, sep = "\t", quote = FALSE, na = "NA")
fwrite(high, out_high, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, out_summary, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Unique missense variants:", nrow(miss_unique), "\n")
cat("High priority missense variants:", nrow(high), "\n")
cat("Output unique:", out_unique, "\n")
cat("Output high:", out_high, "\n")
cat("Summary:", out_summary, "\n")
print(summary)
