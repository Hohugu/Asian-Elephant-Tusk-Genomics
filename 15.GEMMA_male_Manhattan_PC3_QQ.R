#!/usr/bin/env Rscript

library(data.table)

assoc_file <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/result/GEMMA_male/output/tusk_GEMMA_LMM_male_PC3.assoc.txt"

fig_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/figures_GEMMA_male"
tab_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/tables_GEMMA_male"

dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(tab_dir, showWarnings=FALSE, recursive=TRUE)

cat("Loading data...\n")

gwas <- fread(
  assoc_file,
  select=c("chr","rs","ps","p_wald")
)

gwas <- gwas[
  !is.na(p_wald) &
  is.finite(p_wald) &
  p_wald > 0 &
  p_wald <= 1
]

# Main chromosomes only
gwas <- gwas[grepl("^CM044", chr)]

gwas[, chr_num := as.integer(sub("^CM044([0-9]+)\\.1$", "\\1", chr))]
gwas <- gwas[chr_num >= 20 & chr_num <= 48]

# Rename CM044020-CM044048 as chromosomes 1-29
gwas[, CHR := chr_num - 19]
gwas[, logp := -log10(p_wald)]

# Normalize position within each chromosome so all chr have equal width
gwas[, chr_max := max(ps), by=CHR]
gwas[, x := CHR + (ps / chr_max - 0.5) * 0.8]

bonf <- 0.05 / nrow(gwas)
bonf_log <- -log10(bonf)

sig <- gwas[p_wald < bonf]
top10 <- sig[order(p_wald)][1:min(10, .N)]

fwrite(
  sig,
  file.path(tab_dir, "GEMMA_male_PC3_Bonferroni_publication.tsv"),
  sep="\t"
)

cat("SNPs:", nrow(gwas), "\n")
cat("Bonferroni:", bonf, "\n")
cat("-log10 Bonferroni:", bonf_log, "\n")
cat("Significant SNPs:", nrow(sig), "\n")

png(
  file.path(fig_dir, "GEMMA_male_PC3_Manhattan_publication.png"),
  width=3200,
  height=1500,
  res=200
)

par(mar=c(5,5,4,2))

plot(
  gwas$x,
  gwas$logp,
  pch=20,
  cex=0.12,
  col=ifelse(gwas$CHR %% 2 == 0, "grey75", "grey40"),
  xaxt="n",
  xlim=c(0.5,29.5),
  ylim=c(0,18),
  xlab="Chromosome",
  ylab=expression(-log[10](italic(P)[wald])),
  main="GEMMA male-only GWAS PC3"
)

points(
  sig$x,
  sig$logp,
  pch=20,
  cex=0.7,
  col="red"
)

abline(
  h=bonf_log,
  col="red",
  lty=2,
  lwd=2
)

axis(
  1,
  at=1:29,
  labels=1:29,
  las=1,
  cex.axis=0.75
)

axis(
  2,
  at=seq(0,18,2),
  las=1
)

text(
  top10$x,
  top10$logp + 0.35,
  labels=paste0(top10$chr, ":", top10$ps),
  cex=0.42,
  srt=35,
  adj=0
)

legend(
  "topright",
  legend=c("All SNPs", "Bonferroni significant"),
  pch=20,
  col=c("grey55","red"),
  bty="n"
)

dev.off()

# Zoom version
zoom <- gwas[logp > 6]
top15 <- sig[order(p_wald)][1:min(15, .N)]

png(
  file.path(fig_dir, "GEMMA_male_PC3_Manhattan_publication_zoom.png"),
  width=3200,
  height=1500,
  res=200
)

par(mar=c(5,5,4,2))

plot(
  zoom$x,
  zoom$logp,
  pch=20,
  cex=0.6,
  col=ifelse(zoom$p_wald < bonf, "red", "grey45"),
  xaxt="n",
  xlim=c(0.5,29.5),
  ylim=c(6,18),
  xlab="Chromosome",
  ylab=expression(-log[10](italic(P)[wald])),
  main="GEMMA male-only GWAS PC3: strongest signals"
)

abline(h=bonf_log, col="red", lty=2, lwd=2)

axis(1, at=1:29, labels=1:29, las=1, cex.axis=0.75)
axis(2, at=seq(6,18,2), las=1)

text(
  top15$x,
  top15$logp + 0.3,
  labels=paste0(top15$chr, ":", top15$ps),
  cex=0.42,
  srt=35,
  adj=0
)

dev.off()

cat("Done\n")


