#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

conv_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/final_candidate_regions_summary.tsv"
)

camp_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/Campbell_candidates_selection_signals.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/convergent_regions_distance_to_Campbell_genes.tsv"
)

conv <- fread(conv_file)
camp <- fread(camp_file)

conv <- conv[, .(
  Region_ID,
  Region_chr = CHROM,
  Region_start = START,
  Region_end = END,
  Region_gene = Gene_description
)]

camp <- camp[, .(
  Campbell_name,
  Campbell_gene_id = gene_id,
  Campbell_description = description,
  Campbell_chr = CM_chr,
  Campbell_start = start,
  Campbell_end = end
)]

res <- CJ(
  Region_ID = conv$Region_ID,
  Campbell_gene_id = camp$Campbell_gene_id,
  unique = TRUE
)

res <- merge(res, conv, by = "Region_ID")
res <- merge(res, camp, by = "Campbell_gene_id")

res[, Same_chr := Region_chr == Campbell_chr]

res[, Distance_bp := fifelse(
  !Same_chr,
  NA_real_,
  pmax(Campbell_start - Region_end, Region_start - Campbell_end, 0)
)]

res[, Distance_Mb := Distance_bp / 1e6]

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
        fifelse(Distance_bp <= 5000000, "<=5Mb", ">5Mb")
      )
    )
  )
)]

setorder(res, Region_ID, Same_chr, Distance_bp)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(res[Same_chr == TRUE])
