#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

top_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/priority_regions_deltaAF_top100_variants.tsv"
)

genes_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/convergent_regions_annotation_250kb/all_genes.bed"
)

map_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/CM_to_NC.map"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/priority_regions_deltaAF_top100_annotated.tsv"
)

top <- fread(top_file)
genes <- fread(genes_file, header = FALSE)
map <- fread(map_file, header = FALSE)

setnames(map, c("CM_chr", "NC_chr"))

setnames(
  genes,
  c(
    "NC_chr",
    "Gene_start",
    "Gene_end",
    "Gene_ID",
    "Gene_symbol",
    "Description",
    "Gene_biotype"
  )
)

top <- merge(
  top,
  map,
  by.x = "CHROM",
  by.y = "CM_chr",
  all.x = TRUE
)

top[, SNP_start := POS - 1]
top[, SNP_end := POS]

annot_list <- list()

for (i in seq_len(nrow(top))) {
  snp <- top[i]

  g <- genes[NC_chr == snp$NC_chr]

  overlap <- g[Gene_start <= snp$SNP_end & Gene_end >= snp$SNP_start]

  if (nrow(overlap) > 0) {
    tmp <- cbind(snp, overlap)
    tmp[, Relation := "inside_gene"]
    tmp[, Distance_to_gene_bp := 0]
  } else {
    g[, dist := pmax(Gene_start - snp$SNP_end, snp$SNP_start - Gene_end, 0)]
    nearest <- g[which.min(dist)]
    tmp <- cbind(snp, nearest)
    tmp[, Relation := "nearest_gene"]
    tmp[, Distance_to_gene_bp := dist]
    tmp[, dist := NULL]
  }

  annot_list[[i]] <- tmp
}

res <- rbindlist(annot_list, fill = TRUE)

priority_genes <- c(
  "AMELX", "MEP1A", "MEP1B", "PLA2G7", "PDK3",
  "LOC126077071", "LOC126077054", "LOC126074222"
)

res[, Is_priority_gene := Gene_ID %in% priority_genes | Gene_symbol %in% priority_genes]

setorder(res, -abs_delta_AF, Region_ID, POS)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Input top variants:", nrow(top), "\n")
cat("Annotated rows:", nrow(res), "\n")
cat("Output:", out_file, "\n")

print(res[1:30, .(
  Region_ID, CHROM, POS, SNP,
  AF_TT, AF_TX, delta_AF_TX_minus_TT, abs_delta_AF,
  Gene_ID, Gene_symbol, Description,
  Relation, Distance_to_gene_bp, Is_priority_gene
)])
