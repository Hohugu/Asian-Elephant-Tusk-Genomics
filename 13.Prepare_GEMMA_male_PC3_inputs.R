#!/usr/bin/env Rscript

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fam_f <- file.path(base, "PreGWAS/result/GWAS_qc_named.fam")
pc_f  <- file.path(base, "PreGWAS/result/GWAS_covar_PC3.txt")
out   <- file.path(base, "GWAS/result")

fam <- read.table(fam_f, stringsAsFactors=FALSE)
names(fam) <- c("FID","IID","PID","MID","SEX","PHENO")

cov <- read.table(pc_f, header=TRUE, stringsAsFactors=FALSE)

m <- merge(
  fam[,c("FID","IID")],
  cov,
  by=c("FID","IID"),
  sort=FALSE
)

male <- m[m$sex_bin == 1, ]

gemma_covar <- male[, c(
  "pop_bin",
  "PC1",
  "PC2",
  "PC3"
)]

write.table(
  gemma_covar,
  file.path(out, "GEMMA_male_covar_PC3.txt"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)

cat("Male samples:", nrow(male), "\n")

