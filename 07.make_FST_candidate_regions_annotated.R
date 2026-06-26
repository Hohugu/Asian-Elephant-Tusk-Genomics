library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST/tables"

regions <- fread(
  file.path(base,
            "TT_vs_TX_FST_50kb_candidate_regions_merged.tsv")
)

overlap <- fread(
  file.path(base,
            "FST_regions_gene_overlap.tsv"),
  header = FALSE
)

setnames(overlap,
c("NC_chr_region",
  "Region_start",
  "Region_end",
  "Region_ID",
  "NC_chr_gene",
  "Gene_start",
  "Gene_end",
  "Gene_ID",
  "Gene_symbol",
  "Description")
)

map <- fread(
"/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/tables_GEMMA_male/CM_to_NC.map",
header = FALSE
)

setnames(map, c("CM_chr","NC_chr"))

regions <- merge(
  regions,
  map,
  by.x="CHR",
  by.y="CM_chr",
  all.x=TRUE
)

setnames(regions, "CHR", "CM_chr")

final <- merge(
  regions,
  overlap,
  by="Region_ID",
  all.x=TRUE
)

final[, Chromosome_type := "Autosome"]

final[CM_chr=="CM044047.1",
      Chromosome_type := "X"]

final[CM_chr=="CM044048.1",
      Chromosome_type := "Y"]

setcolorder(
  final,
  c(
    "Region_ID",
    "CM_chr",
    "NC_chr",
    "Chromosome_type",
    "Start",
    "End",
    "Length_bp",
    "N_windows",
    "N_SNPs",
    "Mean_FST",
    "Max_Mean_FST",
    "Max_SNP_FST",
    "Top_SNP",
    "Top_SNP_POS",
    "Gene_ID",
    "Gene_symbol",
    "Description"
  )
)

outfile <- file.path(
  base,
  "TT_vs_TX_FST_50kb_candidate_regions_annotated.tsv"
)

fwrite(
  final,
  outfile,
  sep="\t",
  quote=FALSE,
  na="NA"
)

cat("Regions:", uniqueN(final$Region_ID), "\n")
cat("Annotated rows:", nrow(final), "\n")
cat("Output:", outfile, "\n")
