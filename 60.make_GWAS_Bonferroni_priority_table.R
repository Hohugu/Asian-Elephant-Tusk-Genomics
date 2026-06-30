#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

nearest_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_nearest_TSS_annotation.tsv"
)

cds_effect_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_CDS_SNP_effects.tsv"
)

out_all <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_priority_candidates.tsv"
)

out_top <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_priority_candidates_top.tsv"
)

x <- fread(nearest_file)

optional_cols <- c(
  "Gene_ID",
  "Description",
  "Any_selection_overlap",
  "N_selection_categories",
  "nearest_selection_region",
  "nearest_selection_type",
  "distance_to_nearest_selection_bp"
)

for (cc in optional_cols) {
  if (!cc %in% names(x)) {
    x[, (cc) := NA]
  }
}

x[, selection_flag := FALSE]

if ("Any_selection_overlap" %in% names(x)) {
  x[!is.na(Any_selection_overlap) & Any_selection_overlap == TRUE,
    selection_flag := TRUE]
}

if (file.exists(cds_effect_file)) {

  cds <- fread(cds_effect_file)

  cds_summary <- cds[, .(
    CDS_effect = paste(unique(Predicted_effect), collapse = ";"),
    AA_change = paste(unique(paste0(AA_REF, AA_position, AA_ALT)), collapse = ";"),
    CDS_gene = paste(unique(Gene_ID), collapse = ";"),
    CDS_product = paste(unique(Product), collapse = ";"),
    CDS_REF_match = paste(unique(REF_match), collapse = ";")
  ), by = .(CHROM, POS)]

  x <- merge(
    x,
    cds_summary,
    by = c("CHROM", "POS"),
    all.x = TRUE
  )

} else {

  x[, CDS_effect := NA_character_]
  x[, AA_change := NA_character_]
  x[, CDS_gene := NA_character_]
  x[, CDS_product := NA_character_]
  x[, CDS_REF_match := NA_character_]
}

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
    "cytokinesis",
    "DOCK"
  ),
  collapse = "|"
)

x[, text_for_keyword := paste(
  nearest_gene_TSS,
  nearest_gene_description_TSS,
  Gene_ID,
  Description,
  CDS_gene,
  CDS_product,
  sep = " "
)]

x[, has_biological_keyword := grepl(
  bio_pattern,
  text_for_keyword,
  ignore.case = TRUE
)]

x[, priority_score := 0L]

x[grepl("missense|stop_gained|stop_lost", CDS_effect),
  priority_score := priority_score + 5L]

x[grepl("synonymous", CDS_effect),
  priority_score := priority_score + 1L]

x[regulatory_class == "CDS",
  priority_score := priority_score + 3L]

x[regulatory_class == "exon_non_CDS",
  priority_score := priority_score + 2L]

x[regulatory_class == "promoter_2kb",
  priority_score := priority_score + 4L]

x[regulatory_class == "intragenic_non_exonic",
  priority_score := priority_score + 1L]

x[distance_to_nearest_TSS_bp <= 2000,
  priority_score := priority_score + 3L]

x[distance_to_nearest_TSS_bp > 2000 &
    distance_to_nearest_TSS_bp <= 10000,
  priority_score := priority_score + 2L]

x[distance_to_nearest_TSS_bp > 10000 &
    distance_to_nearest_TSS_bp <= 50000,
  priority_score := priority_score + 1L]

x[selection_flag == TRUE,
  priority_score := priority_score + 3L]

x[GWAS_model == "GEMMA_male_PC3",
  priority_score := priority_score + 2L]

x[has_biological_keyword == TRUE,
  priority_score := priority_score + 2L]

x[, priority_reason := ""]

x[grepl("missense", CDS_effect),
  priority_reason := paste0(priority_reason, "missense_CDS;")]

x[grepl("stop_gained|stop_lost", CDS_effect),
  priority_reason := paste0(priority_reason, "stop_variant;")]

x[grepl("synonymous", CDS_effect),
  priority_reason := paste0(priority_reason, "synonymous_CDS;")]

x[regulatory_class == "promoter_2kb",
  priority_reason := paste0(priority_reason, "promoter_2kb;")]

x[regulatory_class == "exon_non_CDS",
  priority_reason := paste0(priority_reason, "exon_non_CDS;")]

x[regulatory_class == "intragenic_non_exonic",
  priority_reason := paste0(priority_reason, "intragenic_non_exonic;")]

x[distance_to_nearest_TSS_bp <= 2000,
  priority_reason := paste0(priority_reason, "TSS_within_2kb;")]

x[distance_to_nearest_TSS_bp > 2000 &
    distance_to_nearest_TSS_bp <= 10000,
  priority_reason := paste0(priority_reason, "TSS_within_10kb;")]

x[distance_to_nearest_TSS_bp > 10000 &
    distance_to_nearest_TSS_bp <= 50000,
  priority_reason := paste0(priority_reason, "TSS_within_50kb;")]

x[selection_flag == TRUE,
  priority_reason := paste0(priority_reason, "selection_overlap;")]

x[GWAS_model == "GEMMA_male_PC3",
  priority_reason := paste0(priority_reason, "male_model;")]

x[has_biological_keyword == TRUE,
  priority_reason := paste0(priority_reason, "biological_keyword;")]

x[, priority_reason := sub(";$", "", priority_reason)]

x[priority_reason == "",
  priority_reason := "low_priority_by_current_rules"]

keep_cols <- c(
  "GWAS_model",
  "SNP",
  "CHROM",
  "POS",
  "P",
  "regulatory_class",
  "priority_score",
  "priority_reason",
  "CDS_effect",
  "AA_change",
  "CDS_gene",
  "CDS_product",
  "nearest_gene_TSS",
  "nearest_gene_description_TSS",
  "nearest_gene_biotype_TSS",
  "distance_to_nearest_TSS_bp",
  "distance_to_gene_body_bp",
  "SNP_location_relative_to_gene",
  "Any_selection_overlap",
  "N_selection_categories",
  "nearest_selection_region",
  "nearest_selection_type",
  "distance_to_nearest_selection_bp",
  "has_biological_keyword"
)

keep_cols <- keep_cols[keep_cols %in% names(x)]

res <- x[, ..keep_cols]

setorder(
  res,
  -priority_score,
  GWAS_model,
  distance_to_nearest_TSS_bp,
  P
)

top <- res[
  priority_score >= 6 |
    regulatory_class %in% c("CDS", "promoter_2kb") |
    grepl("missense|stop_gained|stop_lost", CDS_effect) |
    Any_selection_overlap == TRUE
]

fwrite(res, out_all, sep = "\t", quote = FALSE, na = "NA")
fwrite(top, out_top, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("All priority table:", out_all, "\n")
cat("Top priority table:", out_top, "\n")
cat("Rows all:", nrow(res), "\n")
cat("Rows top:", nrow(top), "\n")

print(
  res[, .N, by = .(GWAS_model, regulatory_class)][order(GWAS_model, -N)]
)

cat("\nTop priority preview:\n")
print(head(top, 30))
