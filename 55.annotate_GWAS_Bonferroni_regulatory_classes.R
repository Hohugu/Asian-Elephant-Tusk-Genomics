#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

snp_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/ALL_GWAS_Bonferroni_SNPs_annotated_fullGFF.tsv"
)

regdir <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_regulatory_annotation.tsv"
)

summary_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_regulatory_summary.tsv"
)

snps <- fread(snp_file)

snps[, SNP_start := POS - 1]
snps[, SNP_end := POS]

setkey(snps, CHROM, SNP_start, SNP_end)

annotate_overlap <- function(snps, bed_file, label) {
  bed <- fread(bed_file, header = FALSE)
  setnames(bed, c("CHROM", "START", "END"))
  bed[, value := TRUE]
  setkey(bed, CHROM, START, END)

  hits <- foverlaps(
    snps,
    bed,
    by.x = c("CHROM", "SNP_start", "SNP_end"),
    by.y = c("CHROM", "START", "END"),
    nomatch = 0
  )

  unique(hits$SNP)
}

cds_hits <- annotate_overlap(snps, file.path(regdir, "CDS.CM.bed"), "CDS")
exon_hits <- annotate_overlap(snps, file.path(regdir, "exons.CM.bed"), "exon")
gene_hits <- annotate_overlap(snps, file.path(regdir, "genes.CM.bed"), "gene")
promoter_hits <- annotate_overlap(snps, file.path(regdir, "promoters_2kb.CM.bed"), "promoter")

snps[, in_CDS := SNP %in% cds_hits]
snps[, in_exon := SNP %in% exon_hits]
snps[, in_gene := SNP %in% gene_hits]
snps[, in_promoter_2kb := SNP %in% promoter_hits]

snps[, regulatory_class := fifelse(
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

summary <- snps[
  ,
  .N,
  by = .(GWAS_model, regulatory_class)
][order(GWAS_model, -N)]

fwrite(snps, out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(summary)
