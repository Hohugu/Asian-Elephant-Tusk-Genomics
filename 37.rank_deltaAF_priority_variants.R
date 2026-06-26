#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
tab <- file.path(base, "Genetics_Analysis/Functional_annotation/tables")

regions <- c("FST_region_18", "FST_region_9", "FST_region_56")

all <- list()

for (r in regions) {
  tt <- fread(file.path(tab, paste0(r, "_TT_AF.tsv")), header = FALSE)
  tx <- fread(file.path(tab, paste0(r, "_TX_AF.tsv")), header = FALSE)

  setnames(tt, c("CHROM","POS","SNP","REF","ALT","AC_TT","AN_TT","AF_TT"))
  setnames(tx, c("CHROM","POS","SNP","REF","ALT","AC_TX","AN_TX","AF_TX"))

  x <- merge(
    tt,
    tx,
    by = c("CHROM","POS","SNP","REF","ALT"),
    all = TRUE
  )

  x[, Region_ID := r]
  x[, delta_AF_TX_minus_TT := AF_TX - AF_TT]
  x[, abs_delta_AF := abs(delta_AF_TX_minus_TT)]

  all[[r]] <- x
}

res <- rbindlist(all, fill = TRUE)
setorder(res, -abs_delta_AF)

fwrite(
  res,
  file.path(tab, "priority_regions_deltaAF_all_variants.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  res[1:100],
  file.path(tab, "priority_regions_deltaAF_top100_variants.tsv"),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Done\n")
cat("Variants:", nrow(res), "\n")
cat("Output:", file.path(tab, "priority_regions_deltaAF_all_variants.tsv"), "\n")
print(res[1:20, .(
  Region_ID, CHROM, POS, SNP, REF, ALT,
  AF_TT, AF_TX, delta_AF_TX_minus_TT, abs_delta_AF
)])
