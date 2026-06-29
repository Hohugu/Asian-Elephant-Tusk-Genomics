#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

rare_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_differentiated_variants.tsv"
)

priority_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_priority_candidates_top.tsv"
)

selection_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/rawVCF_rare_variants_overlap_selection_GWAS.tsv"
)

regdir <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory"
)

gff_file <- "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

map_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/CM_to_NC.map"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_regulatory_priority_annotation.tsv"
)

summary_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_regulatory_priority_summary.tsv"
)

priority_overlap_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/TT_TX_rawVCF_rare_variants_in_GWAS_priority_genes.tsv"
)

rare <- fread(rare_file)

rare[, VAR_start := POS - 1L]
rare[, VAR_end := POS]

setkey(rare, CHROM, VAR_start, VAR_end)

annotate_overlap <- function(dt, bed_file) {

  bed <- fread(bed_file, header = FALSE)
  setnames(bed, c("CHROM", "START", "END"))
  bed[, value := TRUE]
  setkey(bed, CHROM, START, END)

  hits <- foverlaps(
    dt,
    bed,
    by.x = c("CHROM", "VAR_start", "VAR_end"),
    by.y = c("CHROM", "START", "END"),
    nomatch = 0
  )

  unique(hits$SNP)
}

cds_hits <- annotate_overlap(rare, file.path(regdir, "CDS.CM.bed"))
exon_hits <- annotate_overlap(rare, file.path(regdir, "exons.CM.bed"))
gene_hits <- annotate_overlap(rare, file.path(regdir, "genes.CM.bed"))
promoter_hits <- annotate_overlap(rare, file.path(regdir, "promoters_2kb.CM.bed"))

rare[, in_CDS := SNP %in% cds_hits]
rare[, in_exon := SNP %in% exon_hits]
rare[, in_gene := SNP %in% gene_hits]
rare[, in_promoter_2kb := SNP %in% promoter_hits]

rare[, regulatory_class := fifelse(
  in_CDS, "CDS",
  fifelse(
    in_exon, "exon_non_CDS",
    fifelse(
      in_gene, "intragenic_non_exonic",
      fifelse(
        in_promoter_2kb, "promoter_2kb",
        "intergenic"
      )
    )
  )
)]

if (file.exists(selection_file)) {

  sel <- fread(selection_file)

  keep_sel <- intersect(
    c(
      "SNP",
      "in_convergent_region",
      "convergent_region_id",
      "in_XY_score3_region",
      "XY_score3_region_id",
      "is_GWAS_Bonferroni_SNP"
    ),
    names(sel)
  )

  sel <- sel[, ..keep_sel]

  rare <- merge(
    rare,
    sel,
    by = "SNP",
    all.x = TRUE
  )

} else {

  rare[, in_convergent_region := NA]
  rare[, convergent_region_id := NA_character_]
  rare[, in_XY_score3_region := NA]
  rare[, XY_score3_region_id := NA_character_]
  rare[, is_GWAS_Bonferroni_SNP := NA]
}

for (cc in c("in_convergent_region", "in_XY_score3_region", "is_GWAS_Bonferroni_SNP")) {
  if (cc %in% names(rare)) {
    rare[is.na(get(cc)), (cc) := FALSE]
  }
}

map <- fread(map_file, header = FALSE)
setnames(map, c("CM_chr", "NC_chr"))

gff <- fread(
  cmd = paste("awk '$1 !~ /^#/'", shQuote(gff_file)),
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = ""
)

setnames(
  gff,
  c(
    "NC_chr",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attr"
  )
)

genes <- gff[feature == "gene"]

genes[, Gene_ID := fifelse(
  grepl("Name=", attr),
  sub(".*Name=([^;]+).*", "\\1", attr),
  fifelse(
    grepl("gene=", attr),
    sub(".*gene=([^;]+).*", "\\1", attr),
    NA_character_
  )
)]

genes[, Description := fifelse(
  grepl("description=", attr),
  sub(".*description=([^;]+).*", "\\1", attr),
  Gene_ID
)]

genes[, Description := gsub("%2C", ",", Description)]
genes[, Description := gsub("%20", " ", Description)]

genes[, Gene_biotype := fifelse(
  grepl("gene_biotype=", attr),
  sub(".*gene_biotype=([^;]+).*", "\\1", attr),
  NA_character_
)]

genes <- merge(genes, map, by = "NC_chr", all.x = TRUE)

genes[, CHROM := fifelse(!is.na(CM_chr), CM_chr, NC_chr)]

genes[, TSS := fifelse(strand == "+", start, end)]

genes <- genes[, .(
  CHROM,
  Gene_ID,
  Description,
  Gene_biotype,
  strand,
  gene_start = start,
  gene_end = end,
  TSS
)]

out <- vector("list", nrow(rare))

for (i in seq_len(nrow(rare))) {

  s <- rare[i]
  g <- genes[CHROM == s$CHROM]

  if (nrow(g) == 0) {

    out[[i]] <- cbind(
      s,
      data.table(
        nearest_gene_TSS = NA_character_,
        nearest_gene_description_TSS = NA_character_,
        nearest_gene_biotype_TSS = NA_character_,
        nearest_TSS = NA_integer_,
        distance_to_nearest_TSS_bp = NA_real_,
        distance_to_gene_body_bp = NA_real_,
        SNP_location_relative_to_gene = NA_character_
      )
    )

    next
  }

  g[, distance_to_TSS := abs(TSS - s$POS)]

  g[, distance_to_gene_body := fifelse(
    s$POS >= gene_start & s$POS <= gene_end,
    0L,
    pmin(abs(s$POS - gene_start), abs(s$POS - gene_end))
  )]

  best <- g[which.min(distance_to_TSS)]

  relation <- ifelse(
    s$POS >= best$gene_start & s$POS <= best$gene_end,
    "inside_gene",
    ifelse(
      s$POS < best$gene_start,
      "upstream_of_gene",
      "downstream_of_gene"
    )
  )

  out[[i]] <- cbind(
    s,
    data.table(
      nearest_gene_TSS = best$Gene_ID,
      nearest_gene_description_TSS = best$Description,
      nearest_gene_biotype_TSS = best$Gene_biotype,
      nearest_TSS = best$TSS,
      distance_to_nearest_TSS_bp = as.numeric(best$distance_to_TSS),
      distance_to_gene_body_bp = as.numeric(best$distance_to_gene_body),
      SNP_location_relative_to_gene = relation
    )
  )
}

res <- rbindlist(out, fill = TRUE)

priority <- fread(priority_file)

priority_genes <- unique(na.omit(c(
  priority$nearest_gene_TSS,
  unlist(strsplit(paste(na.omit(priority$CDS_gene), collapse = ";"), ";"))
)))

priority_genes <- priority_genes[
  priority_genes != "" &
    priority_genes != "NA" &
    priority_genes != "<NA>"
]

res[, nearest_gene_is_GWAS_priority := nearest_gene_TSS %in% priority_genes]

res[, rare_variant_priority_reason := ""]

res[regulatory_class == "CDS",
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "CDS;")]

res[regulatory_class == "promoter_2kb",
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "promoter_2kb;")]

res[distance_to_nearest_TSS_bp <= 2000,
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "TSS_within_2kb;")]

res[distance_to_nearest_TSS_bp > 2000 &
    distance_to_nearest_TSS_bp <= 10000,
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "TSS_within_10kb;")]

res[nearest_gene_is_GWAS_priority == TRUE,
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "GWAS_priority_gene;")]

res[in_convergent_region == TRUE,
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "convergent_selection;")]

res[in_XY_score3_region == TRUE,
  rare_variant_priority_reason := paste0(rare_variant_priority_reason, "XY_selection;")]

res[, rare_variant_priority_reason := sub(";$", "", rare_variant_priority_reason)]

res[rare_variant_priority_reason == "",
  rare_variant_priority_reason := "low_priority_by_current_rules"]

summary <- rbindlist(
  list(
    data.table(summary_type = "total", N = nrow(res)),
    data.table(summary_type = "CDS", N = sum(res$regulatory_class == "CDS", na.rm = TRUE)),
    data.table(summary_type = "exon_non_CDS", N = sum(res$regulatory_class == "exon_non_CDS", na.rm = TRUE)),
    data.table(summary_type = "intragenic_non_exonic", N = sum(res$regulatory_class == "intragenic_non_exonic", na.rm = TRUE)),
    data.table(summary_type = "promoter_2kb", N = sum(res$regulatory_class == "promoter_2kb", na.rm = TRUE)),
    data.table(summary_type = "intergenic", N = sum(res$regulatory_class == "intergenic", na.rm = TRUE)),
    data.table(summary_type = "within_2kb_TSS", N = sum(res$distance_to_nearest_TSS_bp <= 2000, na.rm = TRUE)),
    data.table(summary_type = "within_10kb_TSS", N = sum(res$distance_to_nearest_TSS_bp <= 10000, na.rm = TRUE)),
    data.table(summary_type = "nearest_GWAS_priority_gene", N = sum(res$nearest_gene_is_GWAS_priority, na.rm = TRUE)),
    data.table(summary_type = "in_convergent_region", N = sum(res$in_convergent_region == TRUE, na.rm = TRUE)),
    data.table(summary_type = "in_XY_score3_region", N = sum(res$in_XY_score3_region == TRUE, na.rm = TRUE))
  ),
  fill = TRUE
)

setorder(res, -abs_delta_AF)

priority_overlap <- res[
  nearest_gene_is_GWAS_priority == TRUE |
    regulatory_class %in% c("CDS", "promoter_2kb") |
    distance_to_nearest_TSS_bp <= 2000 |
    in_convergent_region == TRUE |
    in_XY_score3_region == TRUE
]

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(priority_overlap, priority_overlap_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Annotated rare variants:", nrow(res), "\n")
cat("Priority rare variants:", nrow(priority_overlap), "\n")
cat("Output:", out_file, "\n")
cat("Priority overlap:", priority_overlap_file, "\n")
print(summary)
