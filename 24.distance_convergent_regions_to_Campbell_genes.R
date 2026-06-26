#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

regions_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD.tsv"
)

campbell_file <- file.path(
  base,
  "GWAS/tables_GEMMA_PC3/Campbell_candidates_v2/GEMMA_PC3_all_Campbell_best_pvalues.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/convergent_regions_distance_to_Campbell_candidates.tsv"
)

regions <- fread(regions_file)
campbell <- fread(campbell_file)

campbell_genes <- unique(
  campbell[
    ,
    .(
      Campbell_gene = Gene,
      Campbell_chr = Gene_chr,
      Campbell_start = Gene_start,
      Campbell_end = Gene_end,
      Campbell_description = Gene_description
    )
  ]
)

res <- CJ(
  Region_ID = regions$Region_ID,
  Campbell_gene = campbell_genes$Campbell_gene,
  unique = TRUE
)

res <- merge(res, regions, by = "Region_ID")
res <- merge(res, campbell_genes, by = "Campbell_gene")

res[, Distance_bp := fifelse(
  CHROM != Campbell_chr,
  NA_real_,
  pmax(Campbell_start - END, START - Campbell_end, 0)
)]

res[, Same_chr := CHROM == Campbell_chr]

res[, Distance_class := fifelse(
  !Same_chr,
  "different_chromosome",
  fifelse(
    Distance_bp == 0,
    "overlap",
    fifelse(
      Distance_bp <= 500000,
      "<=500kb",
      fifelse(
        Distance_bp <= 1000000,
        "<=1Mb",
        fifelse(
          Distance_bp <= 2000000,
          "<=2Mb",
          ">2Mb"
        )
      )
    )
  )
)]

setorder(res, Same_chr, Distance_bp)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(res[Same_chr == TRUE, .(Region_ID, CHROM, START, END, Campbell_gene, Campbell_start, Campbell_end, Distance_bp, Distance_class)])
