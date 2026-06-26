#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

campbell_gff <- file.path(base, "GWAS/tables_GEMMA_PC3/Campbell_candidates/Campbell_candidate_gene_positions.gff")
map_file <- file.path(base, "GWAS/tables_GEMMA_male/CM_to_NC.map")

fst_file <- file.path(base, "Genetics_Analysis/Selection/FST/tables/TT_vs_TX_FST_50kb_windows.tsv")
pi_file  <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_pi_50kb_comparison.tsv")
td_file  <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_TajimaD_50kb_comparison.tsv")
het_file <- file.path(base, "Genetics_Analysis/Diversity/TT_TX_windowed_observed_heterozygosity_50kb.tsv")
ld_file  <- file.path(base, "Genetics_Analysis/LD/tables/microLD_summary_TT_TX.tsv")

out_file <- file.path(base, "Genetics_Analysis/Functional_annotation/tables/Campbell_candidates_selection_signals.tsv")

gff <- fread(campbell_gff, header = FALSE)
setnames(gff, c("NC_chr","source","feature","start","end","score","strand","phase","attr"))

gff[, gene_id := sub(".*Name=([^;]+).*", "\\1", attr)]
gff[, description := sub(".*description=([^;]+).*", "\\1", attr)]
gff[, description := gsub("%2C", ",", description)]

map <- fread(map_file, header = FALSE)
setnames(map, c("CM_chr","NC_chr"))

gff <- merge(gff, map, by = "NC_chr", all.x = TRUE)

fst <- fread(fst_file)
pi  <- fread(pi_file)
td  <- fread(td_file)
het <- fread(het_file)
ld  <- fread(ld_file)

get_window_signal <- function(chr, pos1, pos2) {
  center <- floor((pos1 + pos2) / 2)

  f <- fst[CHR == chr & window_start <= center & window_end >= center]
  p <- pi[CHROM == chr & BIN_START <= center & BIN_END >= center]
  t <- td[CHROM == chr & BIN_START <= center & BIN_END >= center]
  h <- het[CHROM == chr & BIN_START <= center & BIN_END >= center]

  data.table(
    Mean_FST = ifelse(nrow(f) > 0, f$Mean_FST[1], NA),
    Max_FST = ifelse(nrow(f) > 0, f$Max_FST[1], NA),
    delta_pi = ifelse(nrow(p) > 0, p$delta_pi[1], NA),
    TajimaD_TT = ifelse(nrow(t) > 0, t$TajimaD_TT[1], NA),
    TajimaD_TX = ifelse(nrow(t) > 0, t$TajimaD_TX[1], NA),
    delta_TajimaD = ifelse(nrow(t) > 0, t$delta_TajimaD[1], NA),
    Ho_TT = ifelse(nrow(h) > 0, h$Ho_TT[1], NA),
    Ho_TX = ifelse(nrow(h) > 0, h$Ho_TX[1], NA),
    delta_Ho_TX_minus_TT = ifelse(nrow(h) > 0, h$delta_Ho_TX_minus_TT[1], NA)
  )
}

res_list <- list()

for (i in seq_len(nrow(gff))) {
  gene <- gff[i]
  sig <- get_window_signal(gene$CM_chr, gene$start, gene$end)

  res_list[[i]] <- cbind(
    gene[, .(gene_id, description, NC_chr, CM_chr, start, end, strand)],
    sig
  )
}

res <- rbindlist(res_list, fill = TRUE)

res[, Campbell_name := fifelse(grepl("amelogenin", description, ignore.case=TRUE), "AMELX/AMELY",
                        fifelse(grepl("meprin", description, ignore.case=TRUE), "MEP1A",
                        fifelse(grepl("enamelin", description, ignore.case=TRUE), "ENAM",
                        fifelse(grepl("ameloblastin", description, ignore.case=TRUE), "AMBN",
                        fifelse(grepl("amelotin", description, ignore.case=TRUE), "AMTN",
                        fifelse(grepl("odontogenic", description, ignore.case=TRUE), "ODAM", gene_id))))))]

res <- merge(
  res,
  ld[, .(
    Campbell_name = region,
    mean_r2_TT,
    mean_r2_TX,
    delta_mean_r2_TX_minus_TT
  )],
  by = "Campbell_name",
  all.x = TRUE
)

setorder(res, Campbell_name)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(res)
