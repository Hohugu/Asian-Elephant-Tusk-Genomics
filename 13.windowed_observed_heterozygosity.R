#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"
div <- file.path(base, "Genetics_Analysis/Diversity")

vcf <- file.path(div, "GWAS_qc_named.vcf.gz")
tt_samples <- fread(file.path(div, "TT.vcf.samples"), header = FALSE)$V1
tx_samples <- fread(file.path(div, "TX.vcf.samples"), header = FALSE)$V1

out_file <- file.path(div, "TT_TX_windowed_observed_heterozygosity_50kb.tsv")

cmd <- paste("bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'", vcf)

con <- pipe(cmd, "r")
header_samples <- system(paste("bcftools query -l", shQuote(vcf)), intern = TRUE)

tt_idx <- which(header_samples %in% tt_samples)
tx_idx <- which(header_samples %in% tx_samples)

res <- list()
chunk <- 0L

repeat {
  lines <- readLines(con, n = 200000)
  if (length(lines) == 0) break

  dt <- tstrsplit(lines, "\t", fixed = TRUE)
  chrom <- dt[[1]]
  pos <- as.integer(dt[[2]])

  gt <- dt[-c(1,2)]

  get_het <- function(idx) {
   mat <- as.data.table(gt[idx])

   het <- mat[
     ,
     Reduce(
       `+`,
       lapply(.SD, function(x) x %in% c("0/1", "1/0", "0|1", "1|0"))
     )
   ]

   nonmiss <- mat[
     ,
     Reduce(
       `+`,
       lapply(.SD, function(x) !(x %in% c("./.", ".|.", ".")))
     )
   ]

   list(het = het, nonmiss = nonmiss)
 }

  tt <- get_het(tt_idx)
  tx <- get_het(tx_idx)

  tmp <- data.table(
    CHROM = chrom,
    POS = pos,
    BIN_START = floor((pos - 1) / 50000) * 50000 + 1,
    BIN_END = floor((pos - 1) / 50000) * 50000 + 50000,
    TT_het = tt$het,
    TT_nonmiss = tt$nonmiss,
    TX_het = tx$het,
    TX_nonmiss = tx$nonmiss
  )

  win <- tmp[
    ,
    .(
      N_SNPs = .N,
      TT_het = sum(TT_het),
      TT_nonmiss = sum(TT_nonmiss),
      TX_het = sum(TX_het),
      TX_nonmiss = sum(TX_nonmiss)
    ),
    by = .(CHROM, BIN_START, BIN_END)
  ]

  res[[length(res) + 1]] <- win
  chunk <- chunk + 1L
  cat("Processed chunk", chunk, "\n")
}

close(con)

all <- rbindlist(res)

final <- all[
  ,
  .(
    N_SNPs = sum(N_SNPs),
    TT_het = sum(TT_het),
    TT_nonmiss = sum(TT_nonmiss),
    TX_het = sum(TX_het),
    TX_nonmiss = sum(TX_nonmiss)
  ),
  by = .(CHROM, BIN_START, BIN_END)
]

final[, Ho_TT := TT_het / TT_nonmiss]
final[, Ho_TX := TX_het / TX_nonmiss]
final[, delta_Ho_TX_minus_TT := Ho_TX - Ho_TT]

final <- final[N_SNPs >= 10]

fwrite(final, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Windows:", nrow(final), "\n")
cat("Output:", out_file, "\n")
