#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fst_file <- file.path(base, "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_candidate_regions_annotated.tsv")
mixed_file <- file.path(base, "GWAS/tables_GEMMA_PC3/GEMMA_PC3_Bonferroni_genes_summary.tsv")
male_file <- file.path(base, "GWAS/tables_GEMMA_male/GEMMA_male_PC3_significant_genes_summary_README_format.tsv")

outdir <- file.path(base, "Genetics_Analysis/Selection/FST/tables")

fst <- fread(fst_file)
mixed <- fread(mixed_file, header = FALSE, fill = TRUE)
male <- fread(male_file)

mixed <- mixed[, .(
  gene = V1,
  symbol = V2,
  nc_chr = V5,
  cm_chr = V6,
  start = V7,
  end = V8,
  best_snp = V9,
  description = V10,
  best_p = V11
)]

fst_genes <- fst[
  !is.na(Gene_ID) & Gene_ID != "NA",
  .(
    gene = Gene_ID[1],
    symbol_fst = Gene_symbol[1],
    fst_region = paste(unique(Region_ID), collapse = ";"),
    fst_chr = paste(unique(CM_chr), collapse = ";"),
    chr_type = paste(unique(Chromosome_type), collapse = ";"),
    max_mean_fst = max(Max_Mean_FST, na.rm = TRUE),
    max_snp_fst = max(Max_SNP_FST, na.rm = TRUE),
    fst_top_snp = Top_SNP[which.max(Max_SNP_FST)][1],
    description_fst = Description[1]
  ),
  by = Gene_ID
][, Gene_ID := NULL]

mixed_genes <- mixed[
  ,
  .(
    gene = gene,
    symbol_mixed = symbol,
    mixed_chr = cm_chr,
    mixed_start = start,
    mixed_end = end,
    mixed_best_snp = best_snp,
    mixed_description = description,
    mixed_best_p = best_p
  )
]

male_genes <- male[
  ,
  .(
    gene = `Gene name`,
    symbol_male = `Gene symbol`,
    male_chr = `CM chr`,
    male_start = Start,
    male_end = End,
    male_best_snp = `Best SNP`,
    male_description = Description,
    male_best_p = `Best p-value`
  )
]

all_genes <- unique(c(
  fst_genes$gene,
  mixed_genes$gene,
  male_genes$gene
))

res <- data.table(gene = all_genes)

res <- merge(res, fst_genes, by = "gene", all.x = TRUE)
res <- merge(res, mixed_genes, by = "gene", all.x = TRUE)
res <- merge(res, male_genes, by = "gene", all.x = TRUE)

res[, in_fst := !is.na(fst_region)]
res[, in_gwas_mixed := !is.na(mixed_best_snp)]
res[, in_gwas_male := !is.na(male_best_snp)]

res[, evidence := fifelse(
  in_fst & in_gwas_mixed & in_gwas_male, "FST + mixed GWAS + male GWAS",
  fifelse(
    in_fst & in_gwas_mixed, "FST + mixed GWAS",
    fifelse(
      in_fst & in_gwas_male, "FST + male GWAS",
      fifelse(
        in_gwas_mixed & in_gwas_male, "mixed GWAS + male GWAS",
        fifelse(
          in_fst, "FST only",
          fifelse(in_gwas_mixed, "mixed GWAS only", "male GWAS only")
        )
      )
    )
  )
)]

res[, description := fifelse(
  !is.na(description_fst), description_fst,
  fifelse(!is.na(mixed_description), mixed_description, male_description)
)]

setorder(res, -in_fst, -in_gwas_mixed, -in_gwas_male, gene)

outfile <- file.path(outdir, "FST_vs_GWAS_gene_overlap.tsv")
summary_file <- file.path(outdir, "FST_vs_GWAS_gene_overlap_summary.tsv")

fwrite(res, outfile, sep = "\t", quote = FALSE, na = "NA")

summary <- res[, .N, by = evidence][order(-N)]
fwrite(summary, summary_file, sep = "\t", quote = FALSE)

cat("Done\n")
cat("FST genes:", nrow(fst_genes), "\n")
cat("Mixed GWAS genes:", nrow(mixed_genes), "\n")
cat("Male GWAS genes:", nrow(male_genes), "\n")
cat("Total unique genes:", nrow(res), "\n")
cat("Output:", outfile, "\n")
cat("Summary:", summary_file, "\n")
print(summary)
