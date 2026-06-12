#!/usr/bin/env Rscript

assoc_file <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/result/GEMMA_PC3/output/tusk_GEMMA_LMM_PC3.assoc.txt"

fig_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/figures_GEMMA_PC3"
tab_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/tables_GEMMA_PC3"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading GEMMA results...\n")

gwas <- read.table(
  assoc_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

gwas <- gwas[!is.na(gwas$p_wald), ]

cat("Rows:", nrow(gwas), "\n")

gwas$logp <- -log10(gwas$p_wald)

# Manhattan

chr_levels <- unique(gwas$chr)

gwas$CHR_INDEX <- match(gwas$chr, chr_levels)

gwas <- gwas[order(gwas$CHR_INDEX, gwas$ps), ]

gwas$BPcum <- seq_len(nrow(gwas))

png(
  file.path(fig_dir, "GEMMA_PC3_Manhattan.png"),
  width = 1800,
  height = 900
)

plot(
  gwas$BPcum,
  gwas$logp,
  pch = 20,
  cex = 0.4,
  col = ifelse(gwas$CHR_INDEX %% 2 == 0,
               "steelblue",
               "grey40"),
  xlab = "Genome",
  ylab = expression(-log[10](p)),
  main = "GEMMA mixed model GWAS (PC1-PC3)"
)

abline(
  h = -log10(0.05 / nrow(gwas)),
  col = "red",
  lwd = 2,
  lty = 2
)

dev.off()

# QQ plot

obs <- sort(gwas$p_wald)
exp <- ppoints(length(obs))

png(
  file.path(fig_dir, "GEMMA_PC3_QQ.png"),
  width = 900,
  height = 900
)

plot(
  -log10(exp),
  -log10(obs),
  pch = 20,
  cex = 0.5,
  xlab = "Expected -log10(P)",
  ylab = "Observed -log10(P)",
  main = "QQ plot GEMMA PC3"
)

abline(0,1,col="red",lwd=2)

dev.off()

# Lambda

chisq <- qchisq(1 - gwas$p_wald, 1)

lambda <- median(chisq, na.rm = TRUE) /
          qchisq(0.5, 1, lower.tail = FALSE)

write.table(
  data.frame(lambda=lambda),
  file.path(tab_dir, "lambda_GEMMA_PC3.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# Top SNPs

top100 <- gwas[order(gwas$p_wald), ][1:100, ]

write.table(
  top100,
  file.path(tab_dir, "Top100_GEMMA_PC3.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("Finished\n")
cat("Lambda =", lambda, "\n")
