#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

conv_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/FST_pi_TajimaD_convergence.tsv"
)

annot_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/convergent_regions_annotation_250kb/convergent_regions_250kb_gene_overlap.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD.tsv"
)

conv <- fread(conv_file)
annot <- fread(annot_file, header = FALSE)

setnames(
  annot,
  c(
    "NC_region_chr", "Region_250kb_start", "Region_250kb_end",
    "Region_ID", "Core_window", "Core_start", "Core_end",
    "Gene_chr", "Gene_start", "Gene_end",
    "Gene_ID", "Gene_symbol", "Description", "Gene_biotype"
  )
)

annot_summary <- annot[
  ,
  .(
    Nearby_genes = paste(unique(Gene_ID), collapse = ";"),
    Nearby_symbols = paste(unique(Gene_symbol), collapse = ";"),
    Nearby_descriptions = paste(unique(Description), collapse = ";"),
    Protein_coding_genes = paste(unique(Gene_ID[Gene_biotype == "protein_coding"]), collapse = ";"),
    Protein_coding_descriptions = paste(unique(Description[Gene_biotype == "protein_coding"]), collapse = ";"),
    N_nearby_genes = uniqueN(Gene_ID),
    N_protein_coding = uniqueN(Gene_ID[Gene_biotype == "protein_coding"])
  ),
  by = Region_ID
]

regions <- unique(
  conv[
    ,
    .(
      Region_ID,
      CHROM,
      START,
      END,
      Chromosome_type,
      Mean_FST,
      Max_Mean_FST,
      Max_SNP_FST,
      Top_SNP,
      Top_SNP_POS,
      delta_pi,
      TajimaD_TT,
      TajimaD_TX,
      delta_TajimaD
    )
  ]
)

res <- merge(regions, annot_summary, by = "Region_ID", all.x = TRUE)

res[, Direction := fifelse(
  delta_pi < 0 & delta_TajimaD < 0,
  "TX shows reduced diversity and lower TajimaD",
  "direction unclear"
)]

res[, Evidence_score := 0L]
res[!is.na(Max_SNP_FST), Evidence_score := Evidence_score + 1L]
res[delta_pi < 0, Evidence_score := Evidence_score + 1L]
res[delta_TajimaD < 0, Evidence_score := Evidence_score + 1L]

setorder(res, -Evidence_score, delta_pi)

fwrite(res, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Integrated candidate regions:", nrow(res), "\n")
cat("Output:", out_file, "\n")
print(res)
