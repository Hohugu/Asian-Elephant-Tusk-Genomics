#!/usr/bin/env Rscript

base_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025"
pheno_f <- file.path(base_dir, "data/PHENOTYPE/GWAS_complete_sample.txt")
pca_f   <- file.path(base_dir, "PreGWAS/result/GWAS_qc_named_pca.eigenvec")
out_dir <- file.path(base_dir, "PreGWAS/result")

pheno <- read.table(
  pheno_f,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

pheno <- pheno[, 1:7]
names(pheno) <- c("IID", "tusk", "tusk_num", "sex", "sex_bin", "population", "pop_bin")

pca <- read.table(pca_f, header = FALSE, stringsAsFactors = FALSE)
names(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

pheno$FID <- pheno$IID

# Phénotype binaire PLINK :
# TT = 2 = cas
# TX = 1 = contrôle
gwas_pheno <- pheno[, c("FID", "IID", "tusk_num")]

# Covariables
covar <- merge(
  pheno[, c("FID", "IID", "sex_bin", "pop_bin")],
  pca[, c("FID", "IID", paste0("PC", 1:10))],
  by = c("FID", "IID")
)

write.table(
  gwas_pheno,
  file.path(out_dir, "GWAS_pheno.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  covar[, c("FID", "IID", "sex_bin", "pop_bin", paste0("PC", 1:5))],
  file.path(out_dir, "GWAS_covar_PC5.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  covar[, c("FID", "IID", "sex_bin", "pop_bin", paste0("PC", 1:10))],
  file.path(out_dir, "GWAS_covar_PC10.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Files written:\n")
cat(file.path(out_dir, "GWAS_pheno.txt"), "\n")
cat(file.path(out_dir, "GWAS_covar_PC5.txt"), "\n")
cat(file.path(out_dir, "GWAS_covar_PC10.txt"), "\n")
