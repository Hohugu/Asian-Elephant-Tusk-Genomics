#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
rvdir <- file.path(base, "Genetics_Analysis/Functional_annotation/rare_variants")

tt_file <- file.path(rvdir, "TT.frq")
tx_file <- file.path(rvdir, "TX.frq")

out_rare <- file.path(rvdir, "TT_TX_rare_differentiated_variants.tsv")
out_summary <- file.path(rvdir, "TT_TX_rare_variant_summary.tsv")

cmd <- paste(
  "paste",
  shQuote(tt_file),
  shQuote(tx_file),
  "| awk 'NR>1 {",
  "chrom=$1; pos=$2;",
  "an_tt=$4; af_tt=$6;",
  "an_tx=$10; af_tx=$12;",
  "afg=((af_tt*an_tt)+(af_tx*an_tx))/(an_tt+an_tx);",
  "maf=(afg<1-afg?afg:1-afg);",
  "delta=af_tx-af_tt;",
  "absd=(delta<0?-delta:delta);",
  "class=\"not_rare\";",
  "if (maf<=0.05 && af_tx>=0.10 && af_tt<=0.02) class=\"TX_enriched_rare\";",
  "else if (maf<=0.05 && af_tt>=0.10 && af_tx<=0.02) class=\"TT_enriched_rare\";",
  "else if (maf<=0.05) class=\"rare_not_strongly_differentiated\";",
  "if (class==\"TX_enriched_rare\" || class==\"TT_enriched_rare\")",
  "print chrom,pos,an_tt,af_tt,an_tx,af_tx,afg,maf,delta,absd,class;",
  "}'"
)

rare <- fread(
  cmd = cmd,
  header = FALSE
)

if (nrow(rare) == 0) {
  rare <- data.table(
    CHROM=character(), POS=integer(), AN_TT=numeric(), AF_TT=numeric(),
    AN_TX=numeric(), AF_TX=numeric(), AF_global=numeric(), MAF_global=numeric(),
    delta_AF_TX_minus_TT=numeric(), abs_delta_AF=numeric(), rare_class=character()
  )
} else {
  setnames(
    rare,
    c("CHROM","POS","AN_TT","AF_TT","AN_TX","AF_TX",
      "AF_global","MAF_global","delta_AF_TX_minus_TT","abs_delta_AF","rare_class")
  )
  rare[, SNP := paste0(CHROM, ":", POS)]
  setcolorder(rare, c("SNP","CHROM","POS","AN_TT","AF_TT","AN_TX","AF_TX",
                      "AF_global","MAF_global","delta_AF_TX_minus_TT",
                      "abs_delta_AF","rare_class"))
  setorder(rare, -abs_delta_AF)
}

summary <- rare[, .N, by = rare_class][order(-N)]

fwrite(rare, out_rare, sep = "\t", quote = FALSE, na = "NA")
fwrite(summary, out_summary, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Rare differentiated:", nrow(rare), "\n")
cat("Output:", out_rare, "\n")
print(summary)
