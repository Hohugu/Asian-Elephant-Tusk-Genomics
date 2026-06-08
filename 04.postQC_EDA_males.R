#!/usr/bin/env Rscript

base_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025"
res_dir  <- file.path(base_dir, "PreGWAS/result")
pheno_f  <- file.path(base_dir, "data/PHENOTYPE/GWAS_only_male.txt")
fig_dir  <- file.path(base_dir, "PreGWAS/figures_postQCmales")
tab_dir  <- file.path(base_dir, "PreGWAS/tables_postQCmales")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

prefix <- file.path(res_dir, "GWAS_male_qc_maf001_geno005_mind005")

pheno <- read.table(
  pheno_f,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  quote = ""
)

pheno <- pheno[, c(1,2,3,6,7)]
names(pheno) <- c("IID", "tusk", "tusk_num", "population", "pop_bin")
pheno$IID <- as.character(pheno$IID)

freq  <- read.table(
  paste0(prefix, ".freq.frq"),
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

imiss <- read.table(
  paste0(prefix, ".missing.imiss"),
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

lmiss <- read.table(
  paste0(prefix, ".missing.lmiss"),
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

het   <- read.table(
  paste0(prefix, ".het.het"),
  header = TRUE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

bim   <- read.table(
  paste0(prefix, ".bim"),
  header = FALSE,
  fill = TRUE,
  stringsAsFactors = FALSE
)

names(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
imiss$IID <- as.character(imiss$IID)
het$IID <- as.character(het$IID)

sample_qc <- merge(imiss, pheno, by = "IID", all.x = TRUE)
het_qc <- merge(het, pheno, by = "IID", all.x = TRUE)

summary_pheno <- as.data.frame(table(pheno$tusk, pheno$population))
names(summary_pheno) <- c("tusk", "population", "N")

summary_variants <- data.frame(
  n_variants = nrow(freq),
  mean_MAF = mean(freq$MAF, na.rm = TRUE),
  median_MAF = median(freq$MAF, na.rm = TRUE),
  mean_variant_missing = mean(lmiss$F_MISS, na.rm = TRUE),
  median_variant_missing = median(lmiss$F_MISS, na.rm = TRUE),
  mean_sample_missing = mean(imiss$F_MISS, na.rm = TRUE),
  median_sample_missing = median(imiss$F_MISS, na.rm = TRUE),
  mean_F = mean(het$F, na.rm = TRUE),
  median_F = median(het$F, na.rm = TRUE)
)

snps_per_contig <- as.data.frame(table(bim$CHR))
names(snps_per_contig) <- c("CHR", "N_SNP")
snps_per_contig <- snps_per_contig[order(snps_per_contig$N_SNP, decreasing = TRUE), ]

write.table(summary_pheno, file.path(tab_dir, "phenotype_male_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(sample_qc[order(sample_qc$F_MISS, decreasing = TRUE), ],
            file.path(tab_dir, "sample_missingness_with_pheno_male.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(summary_variants, file.path(tab_dir, "variant_global_summary_male.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(snps_per_contig, file.path(tab_dir, "snps_per_contig_male.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(het_qc[order(het_qc$F), ],
            file.path(tab_dir, "heterozygosity_with_pheno_male.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

png(file.path(fig_dir, "01_tusk_counts.png"), width = 900, height = 700)
barplot(table(pheno$tusk),
        main = "Tusk phenotype counts",
        xlab = "Phenotype",
        ylab = "N individuals")
dev.off()

#png(file.path(fig_dir, "02_tusk_by_sex.png"), width = 900, height = 700)
#barplot(table(pheno$tusk, pheno$sex),
#        beside = TRUE,
#        legend = TRUE,
#        main = "Tusk phenotype by sex",
#        xlab = "Sex",
#        ylab = "N individuals")
#dev.off()

png(file.path(fig_dir, "03_tusk_by_population.png"), width = 900, height = 700)
barplot(table(pheno$tusk, pheno$population),
        beside = TRUE,
        legend = TRUE,
        main = "Tusk phenotype by population",
        xlab = "Population",
        ylab = "N individuals")
dev.off()

png(file.path(fig_dir, "04_maf_distribution.png"), width = 900, height = 700)
hist(freq$MAF,
     breaks = 100,
     main = "MAF distribution",
     xlab = "Minor allele frequency",
     ylab = "N variants")
dev.off()

png(file.path(fig_dir, "05_sample_missingness.png"), width = 900, height = 700)
hist(imiss$F_MISS,
     breaks = 50,
     main = "Sample missingness",
     xlab = "F_MISS",
     ylab = "N individuals")
dev.off()

png(file.path(fig_dir, "06_variant_missingness.png"), width = 900, height = 700)
hist(lmiss$F_MISS,
     breaks = 100,
     main = "Variant missingness",
     xlab = "F_MISS",
     ylab = "N variants")
dev.off()

png(file.path(fig_dir, "07_heterozygosity_F.png"), width = 900, height = 700)
hist(het$F,
     breaks = 50,
     main = "Heterozygosity / inbreeding coefficient",
     xlab = "F",
     ylab = "N individuals")
dev.off()

png(file.path(fig_dir, "08_heterozygosity_by_population.png"), width = 900, height = 700)
boxplot(F ~ population,
        data = het_qc,
        main = "Heterozygosity by population",
        xlab = "Population",
        ylab = "F")
dev.off()

png(file.path(fig_dir, "09_heterozygosity_by_tusk.png"), width = 900, height = 700)
boxplot(F ~ tusk,
        data = het_qc,
        main = "Heterozygosity by tusk",
        xlab = "Tusk Type",
        ylab = "F")
dev.off()

#png(file.path(fig_dir, "10_hwe_pvalue_distribution.png"), width = 900, height = 700)
#hist(hwe$P,
#     breaks = 100,
#     main = "Hardy-Weinberg p-value distribution",
#     xlab = "HWE P-value",
#     ylab = "N variants")
#dev.off()

top30 <- head(snps_per_contig, 30)

png(file.path(fig_dir, "11_snps_per_contig_top30.png"), width = 1000, height = 900)
barplot(top30$N_SNP,
        names.arg = top30$CHR,
        las = 2,
        main = "Top 30 contigs by SNP count",
        xlab = "Contig",
        ylab = "N SNPs",
        cex.names = 0.6)
dev.off()

cat("EDA/QC terminé.\n")
cat("Figures :", fig_dir, "\n")
cat("Tables  :", tab_dir, "\n")
