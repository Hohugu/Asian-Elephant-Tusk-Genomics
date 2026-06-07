#!/usr/bin/env Rscript

base_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025"
res_dir  <- file.path(base_dir, "PreGWAS/result")
pheno_f  <- file.path(base_dir, "data/PHENOTYPE/GWAS_complete_sample.txt")
fig_dir  <- file.path(base_dir, "PreGWAS/figures_PCA")
tab_dir  <- file.path(base_dir, "PreGWAS/tables_PCA")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

pca_f <- file.path(res_dir, "GWAS_qc_named_pca.eigenvec")
eig_f <- file.path(res_dir, "GWAS_qc_named_pca.eigenval")

pca <- read.table(pca_f, header = FALSE, stringsAsFactors = FALSE)
eig <- scan(eig_f)

names(pca) <- c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2)))

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

pca$IID <- as.character(pca$IID)
pheno$IID <- as.character(pheno$IID)

dat <- merge(pca, pheno, by = "IID", all.x = TRUE)

var_exp <- eig / sum(eig) * 100
var_table <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  eigenvalue = eig,
  variance_percent = var_exp
)

write.table(dat, file.path(tab_dir, "PCA_with_pheno.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(var_table, file.path(tab_dir, "PCA_variance_explained.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

plot_pca <- function(color_group, filename, main_title) {
  groups <- as.factor(dat[[color_group]])
  cols <- as.numeric(groups)

  png(file.path(fig_dir, filename), width = 1000, height = 800)

  plot(dat$PC1, dat$PC2,
       col = cols,
       pch = 19,
       xlab = paste0("PC1 (", round(var_exp[1], 2), "%)"),
       ylab = paste0("PC2 (", round(var_exp[2], 2), "%)"),
       main = main_title)

  legend("topright",
         legend = levels(groups),
         col = seq_along(levels(groups)),
         pch = 19,
         bty = "n")

  dev.off()
}

plot_pca("population", "01_PCA_PC1_PC2_population.png", "PCA: PC1 vs PC2 by population")
plot_pca("sex",        "02_PCA_PC1_PC2_sex.png",        "PCA: PC1 vs PC2 by sex")
plot_pca("tusk",       "03_PCA_PC1_PC2_tusk.png",       "PCA: PC1 vs PC2 by tusk")

png(file.path(fig_dir, "04_PCA_variance_explained.png"), width = 1000, height = 800)
barplot(var_exp[1:20],
        names.arg = paste0("PC", 1:20),
        las = 2,
        main = "PCA variance explained",
        xlab = "Principal components",
        ylab = "Variance explained (%)")
dev.off()

cat("PCA plots finished\n")
cat("Figures:", fig_dir, "\n")
cat("Tables :", tab_dir, "\n")
