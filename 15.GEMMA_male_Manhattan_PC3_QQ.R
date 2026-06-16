#!/usr/bin/env Rscript

library(data.table)

assoc_file <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/result/GEMMA_male/output/tusk_GEMMA_LMM_male_PC3.assoc.txt"

fig_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/figures_GEMMA_male"
tab_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/tables_GEMMA_male"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading selected columns...\n")

gwas <- fread(
  assoc_file,
  select = c("chr", "rs", "ps", "p_wald")
)

gwas <- gwas[!is.na(p_wald) & p_wald > 0]

gwas[, logp := -log10(p_wald)]

chr_levels <- unique(gwas$chr)
gwas[, CHR_INDEX := match(chr, chr_levels)]
setorder(gwas, CHR_INDEX, ps)
gwas[, BPcum := .I]

cat("Rows:", nrow(gwas), "\n")

png(file.path(fig_dir, "GEMMA_male_PC3_Manhattan_pwald.png"), width=1800, height=900)

plot(
  gwas$BPcum,
  gwas$logp,
  pch=20,
  cex=0.35,
  col=ifelse(gwas$CHR_INDEX %% 2 == 0, "steelblue", "grey40"),
  xlab="Genome",
  ylab=expression(-log[10](p)),
  main="GEMMA male-only mixed-model GWAS (p_wald)"
)

abline(h=-log10(0.05 / nrow(gwas)), col="red", lwd=2, lty=2)

dev.off()

obs <- sort(gwas$p_wald)
exp <- ppoints(length(obs))

png(file.path(fig_dir, "GEMMA_male_PC3_QQ_pwald.png"), width=900, height=900)

plot(
  -log10(exp),
  -log10(obs),
  pch=20,
  cex=0.45,
  xlab="Expected -log10(P)",
  ylab="Observed -log10(P)",
  main="QQ plot GEMMA male-only"
)

abline(0, 1, col="red", lwd=2)

dev.off()

chisq <- qchisq(1 - gwas$p_wald, 1)
lambda <- median(chisq, na.rm=TRUE) / qchisq(0.5, 1, lower.tail=FALSE)

fwrite(data.table(lambda=lambda), file.path(tab_dir, "lambda_GEMMA_male_PC3.tsv"), sep="\t")

top100 <- gwas[order(p_wald)][1:100]
fwrite(top100, file.path(tab_dir, "Top100_GEMMA_male_PC3.tsv"), sep="\t")

cat("Finished.\n")
cat("Lambda =", lambda, "\n")

