#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fst_file <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_candidate_regions_merged.tsv"
)

mixed_file <- file.path(
  base,
  "GWAS/tables_GEMMA_PC3/GEMMA_PC3_Bonferroni_genes_summary.tsv"
)

male_file <- file.path(
  base,
  "GWAS/tables_GEMMA_male/GEMMA_male_PC3_significant_genes_summary_README_format.tsv"
)

outdir <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables"
)

cat("Loading data\n")

fst <- fread(fst_file)

mixed <- fread(
  mixed_file,
  header = FALSE
)

setnames(
  mixed,
  c(
    "Gene_ID",
    "Gene_symbol",
    "NC_chr",
    "CM_chr",
    "Start",
    "End",
    "Best_SNP",
    "Description",
    "P"
  )
)

male <- fread(male_file)

gwas <- rbindlist(
  list(
    mixed[, .(
      Dataset="GWAS_mixed",
      Gene_ID,
      Gene_symbol,
      CM_chr,
      Start,
      End,
      Best_SNP,
      P
    )],

    male[, .(
      Dataset="GWAS_male",
      Gene_ID=`Gene name`,
      Gene_symbol=`Gene symbol`,
      CM_chr=`CM chr`,
      Start,
      End,
      Best_SNP=`Best SNP`,
      P=`Best p-value`
    )]
  ),
  fill=TRUE
)

results <- list()

for(i in seq_len(nrow(gwas))) {

  g <- gwas[i]

  chr <- g$CM_chr

  reg <- fst[CHR == chr]

  if(nrow(reg)==0)
    next

  for(j in seq_len(nrow(reg))) {

    r <- reg[j]

    distance <- max(
      r$Start - g$End,
      g$Start - r$End,
      0
    )

    results[[length(results)+1]] <- data.table(
      Dataset = g$Dataset,
      Gene_ID = g$Gene_ID,
      Gene_symbol = g$Gene_symbol,
      GWAS_chr = chr,
      GWAS_start = g$Start,
      GWAS_end = g$End,
      FST_region = r$Region_ID,
      FST_start = r$Start,
      FST_end = r$End,
      Max_FST = r$Max_SNP_FST,
      Distance_bp = distance
    )
  }
}

res <- rbindlist(results)

nearest <- res[
  ,
  .SD[which.min(Distance_bp)],
  by = .(Dataset,Gene_ID)
]

nearest[
  ,
  Distance_class :=
    fifelse(
      Distance_bp == 0,
      "Overlap",
      fifelse(
        Distance_bp <= 50000,
        "<=50kb",
        fifelse(
          Distance_bp <= 100000,
          "<=100kb",
          fifelse(
            Distance_bp <= 500000,
            "<=500kb",
            ">500kb"
          )
        )
      )
    )
]

outfile <- file.path(
  outdir,
  "FST_GWAS_distance_overlap.tsv"
)

fwrite(nearest,outfile,sep="\t")

summary <- nearest[
  ,
  .N,
  by = .(Dataset,Distance_class)
]

summary_file <- file.path(
  outdir,
  "FST_GWAS_distance_overlap_summary.tsv"
)

fwrite(summary,summary_file,sep="\t")

cat("\nDone\n")
cat("GWAS loci:", nrow(gwas), "\n")
cat("Output:", outfile, "\n\n")

print(summary)
