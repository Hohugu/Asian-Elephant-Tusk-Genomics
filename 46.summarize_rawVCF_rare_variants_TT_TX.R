#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
rvdir <- file.path(base, "Genetics_Analysis/Functional_annotation/rare_variants_rawVCF")

tt_file <- file.path(rvdir, "TT_raw.frq")
tx_file <- file.path(rvdir, "TX_raw.frq")

tt_clean <- file.path(rvdir, "TT_raw.biallelic.clean.tsv")
tx_clean <- file.path(rvdir, "TX_raw.biallelic.clean.tsv")

out_rare <- file.path(rvdir, "TT_TX_rawVCF_rare_differentiated_variants.tsv")
out_summary <- file.path(rvdir, "TT_TX_rawVCF_rare_variant_summary.tsv")

system(sprintf(
  "awk 'BEGIN{OFS=\"\\t\"} NR>1 && $3==2 {print $1,$2,$4,$6}' %s > %s",
  shQuote(tt_file), shQuote(tt_clean)
))

system(sprintf(
  "awk 'BEGIN{OFS=\"\\t\"} NR>1 && $3==2 {print $1,$2,$4,$6}' %s > %s",
  shQuote(tx_file), shQuote(tx_clean)
))

tt <- fread(tt_clean, header = FALSE)
tx <- fread(tx_clean, header = FALSE)

setnames(tt, c("CHROM", "POS", "AN_TT", "AF_TT"))
setnames(tx, c("CHROM", "POS", "AN_TX", "AF_TX"))

x <- merge(tt, tx, by = c("CHROM", "POS"))

x[, AF_global := ((AF_TT * AN_TT) + (AF_TX * AN_TX)) / (AN_TT + AN_TX)]
x[, MAF_global := pmin(AF_global, 1 - AF_global)]
x[, delta_AF_TX_minus_TT := AF_TX - AF_TT]
x[, abs_delta_AF := abs(delta_AF_TX_minus_TT)]

x[, rare_class := fifelse(
  MAF_global <= 0.05 & AF_TX >= 0.10 & AF_TT <= 0.02,
  "TX_enriched_rare",
  fifelse(
    MAF_global <= 0.05 & AF_TT >= 0.10 & AF_TX <= 0.02,
    "TT_enriched_rare",
    "other"
  )
)]

rare <- x[rare_class != "other"]
rare[, SNP := paste0(CHROM, ":", POS)]

setcolorder(rare, c(
  "SNP", "CHROM", "POS", "AN_TT", "AF_TT", "AN_TX", "AF_TX",
  "AF_global", "MAF_global", "delta_AF_TX_minus_TT",
  "abs_delta_AF", "rare_class"
))

setorder(rare, -abs_delta_AF)

summary <- rare[, .N, by = rare_class][order(-N)]

fwrite(rare, out_rare, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, out_summary, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Biallelic merged variants:", nrow(x), "\n")
cat("Rare differentiated:", nrow(rare), "\n")
print(summary)
