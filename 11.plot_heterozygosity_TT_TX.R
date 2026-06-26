#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

het_file <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_heterozygosity.het")
cluster_file <- file.path(base, "Genetics_Analysis/Population_structure/TT_TX_clusters.txt")

outdir <- file.path(base, "Genetics_Analysis/Diversity")
figdir <- file.path(outdir, "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

het <- fread(het_file)
clusters <- fread(cluster_file, header = FALSE)
setnames(clusters, c("FID", "IID", "group"))

dat <- merge(het, clusters, by = c("FID", "IID"))

dat[, Ho := 1 - (`O(HOM)` / `N(NM)`)]

summary <- dat[
  ,
  .(
    N = .N,
    Mean_Ho = mean(Ho, na.rm = TRUE),
    Median_Ho = median(Ho, na.rm = TRUE),
    SD_Ho = sd(Ho, na.rm = TRUE),
    Mean_F = mean(F, na.rm = TRUE),
    Median_F = median(F, na.rm = TRUE)
  ),
  by = group
]

test_Ho <- wilcox.test(Ho ~ group, data = dat)
test_F <- wilcox.test(F ~ group, data = dat)

fwrite(dat, file.path(outdir, "TT_TX_heterozygosity_with_groups.tsv"), sep = "\t")
fwrite(summary, file.path(outdir, "TT_TX_heterozygosity_summary.tsv"), sep = "\t")

fwrite(
  data.table(
    test = c("Ho_Wilcoxon", "F_Wilcoxon"),
    p_value = c(test_Ho$p.value, test_F$p.value)
  ),
  file.path(outdir, "TT_TX_heterozygosity_tests.tsv"),
  sep = "\t"
)

png(file.path(figdir, "TT_TX_observed_heterozygosity_boxplot.png"), width = 900, height = 900)

boxplot(
  Ho ~ group,
  data = dat,
  ylab = "Observed heterozygosity",
  xlab = "Phenotype group",
  main = "Genome-wide observed heterozygosity",
  outline = FALSE
)

stripchart(
  Ho ~ group,
  data = dat,
  vertical = TRUE,
  method = "jitter",
  pch = 20,
  add = TRUE
)

dev.off()

png(file.path(figdir, "TT_TX_inbreeding_coefficient_F_boxplot.png"), width = 900, height = 900)

boxplot(
  F ~ group,
  data = dat,
  ylab = "Inbreeding coefficient F",
  xlab = "Phenotype group",
  main = "Genome-wide inbreeding coefficient",
  outline = FALSE
)

stripchart(
  F ~ group,
  data = dat,
  vertical = TRUE,
  method = "jitter",
  pch = 20,
  add = TRUE
)

dev.off()

cat("Done\n")
cat("Individuals:", nrow(dat), "\n")
cat("Output dir:", outdir, "\n")
print(summary)
cat("Ho Wilcoxon p =", test_Ho$p.value, "\n")
cat("F Wilcoxon p =", test_F$p.value, "\n")
