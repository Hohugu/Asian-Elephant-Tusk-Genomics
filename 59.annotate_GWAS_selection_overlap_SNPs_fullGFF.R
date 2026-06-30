#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

snp_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/GEMMA_Bonferroni_all_and_male_vs_selection_SNP_level.tsv"
)

gff_file <- "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

map_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/CM_to_NC.map"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/GWAS_selection_overlap_SNPs_gene_annotation_fullGFF.tsv"
)

snps <- fread(snp_file)
snps <- snps[Any_selection_overlap == TRUE]

map <- fread(map_file, header = FALSE)
setnames(map, c("CM_chr", "NC_chr"))

tmp <- fread(
  cmd = paste("awk '$1 !~ /^#/'", shQuote(gff_file)),
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = ""
)
gff <- tmp

gff <- gff[V3 == "gene"]
setnames(gff, c("NC_chr","source","feature","start","end","score","strand","phase","attr"))

gff[, Gene_ID := sub(".*Name=([^;]+).*", "\\1", attr)]
gff[, Description := fifelse(
  grepl("description=", attr),
  sub(".*description=([^;]+).*", "\\1", attr),
  Gene_ID
)]
gff[, Description := gsub("%2C", ",", Description)]
gff[, Gene_biotype := fifelse(
  grepl("gene_biotype=", attr),
  sub(".*gene_biotype=([^;]+).*", "\\1", attr),
  NA_character_
)]

gff <- merge(gff, map, by = "NC_chr", all.x = TRUE)

out <- list()

for (i in seq_len(nrow(snps))) {
  s <- snps[i]

  # Use CM if possible, otherwise direct contig name
  g <- gff[CM_chr == s$CHROM | NC_chr == s$CHROM]

  if (nrow(g) == 0) {
    out[[i]] <- data.table(
      SNP = s$SNP,
      CHROM = s$CHROM,
      POS = s$POS,
      GWAS_model = s$GWAS_model,
      P = s$P,
      N_selection_categories = s$N_selection_categories,
      Gene_ID = NA_character_,
      Description = NA_character_,
      Gene_biotype = NA_character_,
      Distance_bp = NA_real_,
      Inside_gene = NA,
      Note = "No gene found on chromosome/contig in GFF or map"
    )
    next
  }

  g[, distance := fifelse(
    s$POS >= start & s$POS <= end,
    0,
    pmin(abs(s$POS - start), abs(s$POS - end))
  )]

  best <- g[which.min(distance)]

  out[[i]] <- data.table(
    SNP = s$SNP,
    CHROM = s$CHROM,
    POS = s$POS,
    GWAS_model = s$GWAS_model,
    P = s$P,
    N_selection_categories = s$N_selection_categories,
    Gene_ID = best$Gene_ID,
    Description = best$Description,
    Gene_biotype = best$Gene_biotype,
    Distance_bp = best$distance,
    Inside_gene = best$distance == 0,
    Note = "Annotated with full GFF"
  )
}

res <- rbindlist(out, fill = TRUE)
setorder(res, -N_selection_categories, P)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Annotated SNPs:", nrow(res), "\n")
cat("Output:", out_file, "\n")
print(res)
