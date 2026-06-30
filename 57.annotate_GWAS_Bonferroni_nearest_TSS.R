#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

gwas_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_regulatory_annotation.tsv"
)

gff_file <- "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

map_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/CM_to_NC.map"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_nearest_TSS_annotation.tsv"
)

summary_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_nearest_TSS_summary.tsv"
)

snps <- fread(gwas_file)

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

out <- vector("list", nrow(snps))

for (i in seq_len(nrow(snps))) {

  s <- snps[i]

  g <- genes[CHROM == s$CHROM]

  if (nrow(g) == 0) {

    out[[i]] <- cbind(
      s,
      data.table(
        nearest_gene_TSS = NA_character_,
        nearest_gene_description_TSS = NA_character_,
        nearest_gene_biotype_TSS = NA_character_,
        nearest_gene_strand_TSS = NA_character_,
        nearest_gene_start_TSS = NA_integer_,
        nearest_gene_end_TSS = NA_integer_,
        nearest_TSS = NA_integer_,
        distance_to_nearest_TSS_bp = NA_integer_,
        distance_to_gene_body_bp = NA_integer_,
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
      nearest_gene_strand_TSS = best$strand,
      nearest_gene_start_TSS = best$gene_start,
      nearest_gene_end_TSS = best$gene_end,
      nearest_TSS = best$TSS,
      distance_to_nearest_TSS_bp = best$distance_to_TSS,
      distance_to_gene_body_bp = best$distance_to_gene_body,
      SNP_location_relative_to_gene = relation
    )
  )
}

res <- rbindlist(out, fill = TRUE)

res[, distance_to_nearest_TSS_bp := as.numeric(distance_to_nearest_TSS_bp)]
res[, distance_to_gene_body_bp := as.numeric(distance_to_gene_body_bp)]

summary <- res[, .(
  N = as.integer(.N),
  median_distance_TSS = as.numeric(
    median(distance_to_nearest_TSS_bp, na.rm = TRUE)
  ),
  N_within_2kb_TSS = as.integer(
    sum(distance_to_nearest_TSS_bp <= 2000, na.rm = TRUE)
  ),
  N_within_10kb_TSS = as.integer(
    sum(distance_to_nearest_TSS_bp <= 10000, na.rm = TRUE)
  ),
  N_within_50kb_TSS = as.integer(
    sum(distance_to_nearest_TSS_bp <= 50000, na.rm = TRUE)
  )
), by = .(GWAS_model, regulatory_class)]

setorder(summary, GWAS_model, regulatory_class)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, summary_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(summary)
