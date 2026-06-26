#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fst <- fread(file.path(base, "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_windows.tsv"))
pi  <- fread(file.path(base, "Genetics_Analysis/Diversity/TT_TX_pi_50kb_comparison.tsv"))
td  <- fread(file.path(base, "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_comparison.tsv"))
het <- fread(file.path(base, "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv"))

outdir <- file.path(base, "Genetics_Analysis/X_chromosome")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

xy_chr <- c("CM044047.1", "CM044048.1")

fst_xy <- fst[CHR %in% xy_chr]
pi_xy  <- pi[CHROM %in% xy_chr]
td_xy  <- td[CHROM %in% xy_chr]
het_xy <- het[CHROM %in% xy_chr]

fwrite(fst_xy, file.path(outdir, "XY_FST_50kb_windows.tsv"), sep = "\t")
fwrite(pi_xy,  file.path(outdir, "XY_pi_50kb_windows.tsv"), sep = "\t")
fwrite(td_xy,  file.path(outdir, "XY_TajimaD_50kb_windows.tsv"), sep = "\t")
fwrite(het_xy, file.path(outdir, "XY_Ho_50kb_windows.tsv"), sep = "\t")

cat("Done\n")
cat("FST XY windows:", nrow(fst_xy), "\n")
cat("pi XY windows:", nrow(pi_xy), "\n")
cat("TajimaD XY windows:", nrow(td_xy), "\n")
cat("Ho XY windows:", nrow(het_xy), "\n")
cat("Output dir:", outdir, "\n")
