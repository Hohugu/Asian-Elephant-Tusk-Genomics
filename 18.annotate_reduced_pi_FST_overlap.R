#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

overlap_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_reduced_pi_overlapping_FST.tsv"
)

fst_annot_file <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_candidate_regions_annotated.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_reduced_pi_overlapping_FST_annotated.tsv"
)

overlap <- fread(overlap_file)
annot <- fread(fst_annot_file)

annot_sub <- annot[
  ,
  .(
    Region_ID,
    CM_chr,
    NC_chr,
    Chromosome_type,
    Start,
    End,
    Length_bp,
    Mean_FST,
    Max_Mean_FST,
    Max_SNP_FST,
    Top_SNP,
    Top_SNP_POS,
    Gene_ID,
    Gene_symbol,
    Description
  )
]

res <- merge(
  overlap,
  annot_sub,
  by = "Region_ID",
  all.x = TRUE,
  allow.cartesian = TRUE
)

setorder(res, CHROM, START, Gene_ID)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Input overlap regions:", uniqueN(overlap$Region_ID), "\n")
cat("Annotated rows:", nrow(res), "\n")
cat("Output:", out_file, "\n")
