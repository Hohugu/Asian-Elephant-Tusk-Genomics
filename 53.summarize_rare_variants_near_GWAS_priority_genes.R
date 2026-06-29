#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

rare_annot_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_regulatory_priority_annotation.tsv"
)

gwas_priority_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_priority_candidates_top.tsv"
)

out_all <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_near_GWAS_priority_genes.tsv"
)

out_high <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_near_GWAS_priority_genes_high_priority.tsv"
)

out_gene_summary <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_near_GWAS_priority_genes_by_gene.tsv"
)

out_summary <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_near_GWAS_priority_genes_summary.tsv"
)

rare <- fread(rare_annot_file)
priority <- fread(gwas_priority_file)

to_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  y <- tolower(as.character(x))
  y %in% c("true", "t", "1", "yes", "y")
}

needed_rare <- c(
  "nearest_gene_is_GWAS_priority",
  "nearest_gene_TSS",
  "distance_to_nearest_TSS_bp",
  "regulatory_class",
  "in_convergent_region",
  "in_XY_score3_region"
)

for (cc in needed_rare) {
  if (!cc %in% names(rare)) {
    rare[, (cc) := NA]
  }
}

rare[, nearest_gene_is_GWAS_priority_bool := to_bool(nearest_gene_is_GWAS_priority)]
rare[, in_convergent_region_bool := to_bool(in_convergent_region)]
rare[, in_XY_score3_region_bool := to_bool(in_XY_score3_region)]

rare_priority <- rare[
  nearest_gene_is_GWAS_priority_bool == TRUE
]

# -------------------------------------------------------------------
# Build GWAS priority gene metadata
# -------------------------------------------------------------------

priority_gene_tables <- list()

if ("nearest_gene_TSS" %in% names(priority)) {
  priority_gene_tables[["nearest"]] <- priority[
    !is.na(nearest_gene_TSS) &
      nearest_gene_TSS != "" &
      nearest_gene_TSS != "NA",
    .(
      priority_gene = nearest_gene_TSS,
      priority_gene_source = "nearest_TSS",
      GWAS_priority_SNP = SNP,
      GWAS_model,
      GWAS_priority_score = priority_score,
      GWAS_priority_reason = priority_reason,
      GWAS_regulatory_class = regulatory_class,
      GWAS_CDS_effect = if ("CDS_effect" %in% names(priority)) CDS_effect else NA_character_,
      GWAS_AA_change = if ("AA_change" %in% names(priority)) AA_change else NA_character_,
      GWAS_gene_description = nearest_gene_description_TSS
    )
  ]
}

if ("CDS_gene" %in% names(priority)) {
  tmp <- priority[
    !is.na(CDS_gene) &
      CDS_gene != "" &
      CDS_gene != "NA",
    .(
      priority_gene_raw = CDS_gene,
      GWAS_priority_SNP = SNP,
      GWAS_model,
      GWAS_priority_score = priority_score,
      GWAS_priority_reason = priority_reason,
      GWAS_regulatory_class = regulatory_class,
      GWAS_CDS_effect = CDS_effect,
      GWAS_AA_change = AA_change,
      GWAS_gene_description = if ("CDS_product" %in% names(priority)) CDS_product else NA_character_
    )
  ]

  if (nrow(tmp) > 0) {
    tmp <- tmp[
      ,
      .(priority_gene = unlist(strsplit(priority_gene_raw, ";", fixed = TRUE))),
      by = .(
        GWAS_priority_SNP,
        GWAS_model,
        GWAS_priority_score,
        GWAS_priority_reason,
        GWAS_regulatory_class,
        GWAS_CDS_effect,
        GWAS_AA_change,
        GWAS_gene_description
      )
    ]

    tmp[, priority_gene_source := "CDS_gene"]

    priority_gene_tables[["cds"]] <- tmp
  }
}

priority_genes_long <- rbindlist(priority_gene_tables, fill = TRUE)

priority_genes_long <- unique(priority_genes_long[
  !is.na(priority_gene) &
    priority_gene != "" &
    priority_gene != "NA" &
    priority_gene != "<NA>"
])

priority_gene_summary <- priority_genes_long[
  ,
  .(
    GWAS_priority_SNPs = paste(sort(unique(GWAS_priority_SNP)), collapse = ";"),
    GWAS_models = paste(sort(unique(GWAS_model)), collapse = ";"),
    GWAS_priority_gene_sources = paste(sort(unique(priority_gene_source)), collapse = ";"),
    best_GWAS_priority_score = suppressWarnings(max(GWAS_priority_score, na.rm = TRUE)),
    GWAS_priority_reasons = paste(sort(unique(GWAS_priority_reason)), collapse = ";"),
    GWAS_regulatory_classes = paste(sort(unique(GWAS_regulatory_class)), collapse = ";"),
    GWAS_CDS_effects = paste(sort(unique(na.omit(GWAS_CDS_effect))), collapse = ";"),
    GWAS_AA_changes = paste(sort(unique(na.omit(GWAS_AA_change))), collapse = ";"),
    GWAS_gene_descriptions = paste(sort(unique(na.omit(GWAS_gene_description))), collapse = ";")
  ),
  by = priority_gene
]

rare_priority <- merge(
  rare_priority,
  priority_gene_summary,
  by.x = "nearest_gene_TSS",
  by.y = "priority_gene",
  all.x = TRUE
)

# -------------------------------------------------------------------
# Rank rare variants near GWAS priority genes
# -------------------------------------------------------------------

rare_priority[, rare_near_GWAS_priority_score := 0L]

rare_priority[regulatory_class == "CDS",
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 4L]

rare_priority[regulatory_class == "promoter_2kb",
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 4L]

rare_priority[regulatory_class == "exon_non_CDS",
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 2L]

rare_priority[regulatory_class == "intragenic_non_exonic",
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 1L]

rare_priority[distance_to_nearest_TSS_bp <= 2000,
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 3L]

rare_priority[distance_to_nearest_TSS_bp > 2000 &
                distance_to_nearest_TSS_bp <= 10000,
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 2L]

rare_priority[distance_to_nearest_TSS_bp > 10000 &
                distance_to_nearest_TSS_bp <= 50000,
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 1L]

rare_priority[in_convergent_region_bool == TRUE,
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 3L]

rare_priority[in_XY_score3_region_bool == TRUE,
  rare_near_GWAS_priority_score := rare_near_GWAS_priority_score + 2L]

rare_priority[, rare_near_GWAS_priority_reason := ""]

rare_priority[regulatory_class == "CDS",
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "CDS;")]

rare_priority[regulatory_class == "promoter_2kb",
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "promoter_2kb;")]

rare_priority[regulatory_class == "exon_non_CDS",
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "exon_non_CDS;")]

rare_priority[regulatory_class == "intragenic_non_exonic",
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "intragenic_non_exonic;")]

rare_priority[distance_to_nearest_TSS_bp <= 2000,
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "TSS_within_2kb;")]

rare_priority[distance_to_nearest_TSS_bp > 2000 &
                distance_to_nearest_TSS_bp <= 10000,
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "TSS_within_10kb;")]

rare_priority[distance_to_nearest_TSS_bp > 10000 &
                distance_to_nearest_TSS_bp <= 50000,
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "TSS_within_50kb;")]

rare_priority[in_convergent_region_bool == TRUE,
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "convergent_selection;")]

rare_priority[in_XY_score3_region_bool == TRUE,
  rare_near_GWAS_priority_reason := paste0(rare_near_GWAS_priority_reason, "XY_selection;")]

rare_priority[, rare_near_GWAS_priority_reason := sub(";$", "", rare_near_GWAS_priority_reason)]

rare_priority[rare_near_GWAS_priority_reason == "",
  rare_near_GWAS_priority_reason := "nearest_GWAS_priority_gene_only"]

setorder(
  rare_priority,
  -rare_near_GWAS_priority_score,
  nearest_gene_TSS,
  distance_to_nearest_TSS_bp,
  -abs_delta_AF
)

high <- rare_priority[
  rare_near_GWAS_priority_score >= 4 |
    regulatory_class %in% c("CDS", "promoter_2kb") |
    distance_to_nearest_TSS_bp <= 10000 |
    in_convergent_region_bool == TRUE |
    in_XY_score3_region_bool == TRUE
]

# -------------------------------------------------------------------
# Gene-level summary
# -------------------------------------------------------------------

gene_summary <- rare_priority[
  ,
  .(
    N_rare_variants = .N,
    max_abs_delta_AF = max(abs_delta_AF, na.rm = TRUE),
    min_distance_to_TSS_bp = min(distance_to_nearest_TSS_bp, na.rm = TRUE),
    N_CDS = sum(regulatory_class == "CDS", na.rm = TRUE),
    N_exon_non_CDS = sum(regulatory_class == "exon_non_CDS", na.rm = TRUE),
    N_promoter_2kb = sum(regulatory_class == "promoter_2kb", na.rm = TRUE),
    N_intragenic_non_exonic = sum(regulatory_class == "intragenic_non_exonic", na.rm = TRUE),
    N_intergenic = sum(regulatory_class == "intergenic", na.rm = TRUE),
    N_within_2kb_TSS = sum(distance_to_nearest_TSS_bp <= 2000, na.rm = TRUE),
    N_within_10kb_TSS = sum(distance_to_nearest_TSS_bp <= 10000, na.rm = TRUE),
    N_XY_score3 = sum(in_XY_score3_region_bool == TRUE, na.rm = TRUE),
    N_convergent = sum(in_convergent_region_bool == TRUE, na.rm = TRUE),
    rare_classes = paste(sort(unique(rare_class)), collapse = ";"),
    rare_variant_SNPs = paste(sort(unique(SNP)), collapse = ";"),
    nearest_gene_description_TSS = paste(sort(unique(na.omit(nearest_gene_description_TSS))), collapse = ";"),
    nearest_gene_biotype_TSS = paste(sort(unique(na.omit(nearest_gene_biotype_TSS))), collapse = ";"),
    GWAS_priority_SNPs = paste(sort(unique(na.omit(GWAS_priority_SNPs))), collapse = ";"),
    GWAS_models = paste(sort(unique(na.omit(GWAS_models))), collapse = ";"),
    GWAS_regulatory_classes = paste(sort(unique(na.omit(GWAS_regulatory_classes))), collapse = ";"),
    GWAS_CDS_effects = paste(sort(unique(na.omit(GWAS_CDS_effects))), collapse = ";"),
    GWAS_AA_changes = paste(sort(unique(na.omit(GWAS_AA_changes))), collapse = ";"),
    best_GWAS_priority_score = suppressWarnings(max(best_GWAS_priority_score, na.rm = TRUE))
  ),
  by = nearest_gene_TSS
]

setorder(
  gene_summary,
  -N_rare_variants,
  min_distance_to_TSS_bp,
  nearest_gene_TSS
)

summary <- data.table(
  metric = c(
    "rare_variants_nearest_GWAS_priority_gene",
    "high_priority_rare_variants_nearest_GWAS_priority_gene",
    "unique_GWAS_priority_genes_with_rare_variant",
    "CDS",
    "promoter_2kb",
    "exon_non_CDS",
    "intragenic_non_exonic",
    "intergenic",
    "within_2kb_TSS",
    "within_10kb_TSS",
    "in_convergent_region",
    "in_XY_score3_region"
  ),
  N = c(
    nrow(rare_priority),
    nrow(high),
    uniqueN(rare_priority$nearest_gene_TSS),
    sum(rare_priority$regulatory_class == "CDS", na.rm = TRUE),
    sum(rare_priority$regulatory_class == "promoter_2kb", na.rm = TRUE),
    sum(rare_priority$regulatory_class == "exon_non_CDS", na.rm = TRUE),
    sum(rare_priority$regulatory_class == "intragenic_non_exonic", na.rm = TRUE),
    sum(rare_priority$regulatory_class == "intergenic", na.rm = TRUE),
    sum(rare_priority$distance_to_nearest_TSS_bp <= 2000, na.rm = TRUE),
    sum(rare_priority$distance_to_nearest_TSS_bp <= 10000, na.rm = TRUE),
    sum(rare_priority$in_convergent_region_bool == TRUE, na.rm = TRUE),
    sum(rare_priority$in_XY_score3_region_bool == TRUE, na.rm = TRUE)
  )
)

fwrite(rare_priority, out_all, sep = "\t", quote = FALSE, na = "NA")
fwrite(high, out_high, sep = "\t", quote = FALSE, na = "NA")
fwrite(gene_summary, out_gene_summary, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, out_summary, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Rare variants near GWAS priority genes:", nrow(rare_priority), "\n")
cat("High priority rare variants:", nrow(high), "\n")
cat("Genes with rare variants:", uniqueN(rare_priority$nearest_gene_TSS), "\n")
cat("Output all:", out_all, "\n")
cat("Output high:", out_high, "\n")
cat("Output gene summary:", out_gene_summary, "\n")
cat("Summary:", out_summary, "\n")
print(summary)

