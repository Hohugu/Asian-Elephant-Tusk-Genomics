#!/usr/bin/env Rscript

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

fam_f <- file.path(base, "PreGWAS/result/GWAS_qc_named.fam")
pc_f  <- file.path(base, "PreGWAS/result/GWAS_covar_PC5.txt")
out   <- file.path(base, "GWAS/result")

fam <- read.table(fam_f, stringsAsFactors = FALSE)
names(fam) <- c("FID","IID","PID","MID","SEX","PHENO")

cov <- read.table(pc_f, header = TRUE, stringsAsFactors = FALSE)

m <- merge(fam[,c("FID","IID")], cov, by=c("FID","IID"), sort=FALSE)

gemma_covar <- m[, c("sex_bin","pop_bin","PC1","PC2","PC3","PC4","PC5")]

write.table(gemma_covar,
            file.path(out, "GEMMA_covar_PC5.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
