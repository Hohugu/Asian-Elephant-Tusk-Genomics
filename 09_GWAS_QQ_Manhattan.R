#!/usr/bin/env Rscript

base_dir <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS"

res_dir <- file.path(base_dir, "result")
fig_dir <- file.path(base_dir, "figures")
tab_dir <- file.path(base_dir, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# LOAD GWAS
# =========================

gwas5 <- read.table(
  file.path(res_dir, "GWAS_logistic_PC5.ADD.clean.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

gwas10 <- read.table(
  file.path(res_dir, "GWAS_logistic_PC10.ADD.clean.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# remove NA p-values just in case
gwas5 <- gwas5[!is.na(gwas5$P), ]
gwas10 <- gwas10[!is.na(gwas10$P), ]

# =========================
# LAMBDA
# =========================

calc_lambda <- function(pvals) {

  chisq <- qchisq(1 - pvals, 1)

  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)

  return(lambda)
}

lambda5 <- calc_lambda(gwas5$P)
lambda10 <- calc_lambda(gwas10$P)

lambda_table <- data.frame(
  Model = c("PC5", "PC10"),
  Lambda = c(lambda5, lambda10)
)

write.table(
  lambda_table,
  file.path(tab_dir, "lambda_inflation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# QQ PLOT
# =========================

make_qq <- function(df, title, outfile) {

  obs <- -log10(sort(df$P))

  exp <- -log10(ppoints(length(obs)))

  png(outfile, width = 1000, height = 1000)

  plot(
    exp,
    obs,
    pch = 20,
    cex = 0.5,
    xlab = "Expected -log10(P)",
    ylab = "Observed -log10(P)",
    main = title
  )

  abline(0, 1, col = "red", lwd = 2)

  dev.off()
}

make_qq(
  gwas5,
  paste0("QQ plot PC5 (lambda=", round(lambda5, 3), ")"),
  file.path(fig_dir, "QQ_PC5.png")
)

make_qq(
  gwas10,
  paste0("QQ plot PC10 (lambda=", round(lambda10, 3), ")"),
  file.path(fig_dir, "QQ_PC10.png")
)

# =========================
# MANHATTAN
# =========================

prepare_manhattan <- function(df) {

  df <- df[order(df$CHR, df$BP), ]

  chr_levels <- unique(df$CHR)

  df$CHR_INDEX <- match(df$CHR, chr_levels)

  chr_sizes <- tapply(df$BP, df$CHR_INDEX, max)

  offsets <- c(0, cumsum(chr_sizes)[-length(chr_sizes)])

  df$BPcum <- df$BP + offsets[df$CHR_INDEX]

  return(df)
}

plot_manhattan <- function(df, title, outfile) {

  df <- prepare_manhattan(df)

  png(outfile, width = 1800, height = 900)

  plot(
    df$BPcum,
    -log10(df$P),
    pch = 20,
    cex = 0.4,
    col = ifelse(df$CHR_INDEX %% 2 == 0, "steelblue", "darkorange"),
    xaxt = "n",
    xlab = "Genome",
    ylab = "-log10(P)",
    main = title
  )

  abline(h = -log10(5e-8), col = "red", lty = 2)

  dev.off()
}

plot_manhattan(
  gwas5,
  "GWAS Manhattan PC5",
  file.path(fig_dir, "Manhattan_PC5.png")
)

plot_manhattan(
  gwas10,
  "GWAS Manhattan PC10",
  file.path(fig_dir, "Manhattan_PC10.png")
)

# =========================
# TOP SNPS
# =========================

top5 <- gwas5[order(gwas5$P), ][1:100, ]
top10 <- gwas10[order(gwas10$P), ][1:100, ]

write.table(
  top5,
  file.path(tab_dir, "Top100_PC5.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  top10,
  file.path(tab_dir, "Top100_PC10.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("GWAS plots finished\n")
cat("Figures:", fig_dir, "\n")
cat("Tables :", tab_dir, "\n")
