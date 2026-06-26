#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

pi_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/TT_TX_pi_50kb_top100_reduced_TX.tsv"
)

fst_file <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_candidate_regions_annotated.tsv"
)

gwas_file <- file.path(
  base,
  "Genetics_Analysis/Selection/FST/tables/FST_GWAS_distance_overlap.tsv"
)

outdir <- file.path(
  base,
  "Genetics_Analysis/Diversity"
)

pi <- fread(pi_file)
fst <- fread(fst_file)

pi[, window :=
      paste(CHROM,BIN_START,BIN_END,sep=":")]

fst_regions <- unique(
  fst[, .(
    CHROM = CM_chr,
    START = Start,
    END = End,
    Region_ID
  )]
)

hits <- list()

for(i in seq_len(nrow(pi))) {

  chr <- pi$CHROM[i]
  s <- pi$BIN_START[i]
  e <- pi$BIN_END[i]

  ov <- fst_regions[
    CHROM==chr &
    START <= e &
    END >= s
  ]

  if(nrow(ov)>0){

    tmp <- copy(ov)

    tmp[, pi_window :=
          paste(chr,s,e,sep=":")]

    tmp[, delta_pi := pi$delta_pi[i]]

    hits[[length(hits)+1]] <- tmp
  }
}

res <- rbindlist(hits, fill=TRUE)

outfile <- file.path(
  outdir,
  "TT_TX_reduced_pi_overlapping_FST.tsv"
)

fwrite(res,outfile,sep="\t")

cat("Reduced-pi windows overlapping FST:",nrow(res),"\n")
cat("Output:",outfile,"\n")
