#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

xy_file <- file.path(
  base,
  "Genetics_Analysis/X_chromosome/XY_selection_signals_top50.tsv"
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
  "Genetics_Analysis/X_chromosome/XY_selection_signals_top50_annotated.tsv"
)

xy <- fread(xy_file)
genes <- fread(genes_file, header = FALSE)
map <- fread(map_file, header = FALSE)

setnames(map, c("CM_chr", "NC_chr"))

setnames(
  genes,
  c("NC_chr", "Gene_start", "Gene_end", "Gene_ID",
    "Gene_symbol", "Description", "Gene_biotype")
)

xy <- merge(
  xy,
  map,
  by.x = "CHROM",
  by.y = "CM_chr",
  all.x = TRUE
)

annot <- list()

for (i in seq_len(nrow(xy))) {
  r <- xy[i]
  g <- genes[NC_chr == r$NC_chr]

  overlap <- g[Gene_start <= r$BIN_END & Gene_end >= r$BIN_START]

  if (nrow(overlap) > 0) {
    tmp <- cbind(r, overlap)
    tmp[, Relation := "overlap"]
    tmp[, Distance_to_gene_bp := 0]
  } else {
    g[, dist := pmax(Gene_start - r$BIN_END, r$BIN_START - Gene_end, 0)]
    nearest <- g[which.min(dist)]
    tmp <- cbind(r, nearest)
    tmp[, Relation := "nearest_gene"]
    tmp[, Distance_to_gene_bp := dist]
    tmp[, dist := NULL]
  }

  annot[[i]] <- tmp
}

res <- rbindlist(annot, fill = TRUE)

setorder(res, -Evidence_score, -Max_FST)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Annotated rows:", nrow(res), "\n")
cat("Output:", out_file, "\n")

print(res[1:30, .(
  Chromosome_type, CHROM, BIN_START, BIN_END,
  Max_FST, delta_pi, delta_Ho_TX_minus_TT,
  Gene_ID, Gene_symbol, Description, Gene_biotype,
  Relation, Distance_to_gene_bp
)])
