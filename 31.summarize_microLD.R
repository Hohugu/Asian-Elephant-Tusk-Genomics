#!/usr/bin/env Rscript

library(data.table)

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
LDDIR <- file.path(BASE, "Genetics_Analysis", "LD", "plink")
OUTDIR <- file.path(BASE, "Genetics_Analysis", "LD", "tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

files <- list.files(LDDIR, pattern = "_LD\\.ld$", full.names = TRUE)

res <- list()

for (f in files) {
  nm <- basename(f)
  nm <- sub("_LD.ld$", "", nm)

  cmd <- sprintf(
    "awk 'NR==1{for(i=1;i<=NF;i++) if($i==\"R2\") col=i; next} {n++; sum+=$col; if($col>=0.8) high++} END{if(n>0) print n, sum/n, high+0; else print 0, \"NA\", 0}' %s",
    shQuote(f)
  )

  x <- strsplit(system(cmd, intern = TRUE), " ")[[1]]

  res[[nm]] <- data.table(
    dataset = nm,
    region = sub("_(TT|TX)$", "", nm),
    group = sub("^.*_(TT|TX)$", "\\1", nm),
    n_pairs = as.numeric(x[1]),
    mean_r2 = as.numeric(x[2]),
    n_r2_over_0.8 = as.numeric(x[3])
  )

  cat("Done:", nm, "\n")
}

out <- rbindlist(res)

wide <- dcast(
  out,
  region ~ group,
  value.var = c("n_pairs", "mean_r2", "n_r2_over_0.8")
)

if ("mean_r2_TX" %in% names(wide) & "mean_r2_TT" %in% names(wide)) {
  wide[, delta_mean_r2_TX_minus_TT := mean_r2_TX - mean_r2_TT]
}

if ("n_r2_over_0.8_TX" %in% names(wide) & "n_r2_over_0.8_TT" %in% names(wide)) {
  wide[, delta_highLD_TX_minus_TT := n_r2_over_0.8_TX - n_r2_over_0.8_TT]
}

setorder(wide, -delta_mean_r2_TX_minus_TT)

fwrite(out, file.path(OUTDIR, "microLD_summary_long.tsv"), sep = "\t")
fwrite(wide, file.path(OUTDIR, "microLD_summary_TT_TX.tsv"), sep = "\t")

cat("Done\n")
cat("Output:", file.path(OUTDIR, "microLD_summary_TT_TX.tsv"), "\n")
print(wide)
