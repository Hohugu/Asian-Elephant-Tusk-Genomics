#!/usr/bin/env Rscript

library(data.table)
options(scipen = 999)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025"

integrated_file <- file.path(
  base,
  "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv"
)

deltaaf_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/priority_regions_deltaAF_top100_annotated.tsv"
)

ld_file <- file.path(
  base,
  "Genetics_Analysis/LD/tables/microLD_summary_TT_TX.tsv"
)

out_file <- file.path(
  base,
  "Genetics_Analysis/Functional_annotation/tables/final_candidate_regions_summary.tsv"
)

x <- fread(integrated_file)
af <- fread(deltaaf_file)

# Remove duplicated column names from annotated deltaAF table
dup <- duplicated(names(af))
if (any(dup)) {
  af <- af[, !dup, with = FALSE]
}

ld <- fread(ld_file)

top_af <- af[
  ,
  .SD[which.max(abs_delta_AF)],
  by = Region_ID
]

ld_sub <- ld[
  ,
  .(
    Region_ID = region,
    mean_r2_TT,
    mean_r2_TX,
    delta_mean_r2_TX_minus_TT,
    n_r2_over_0.8_TT,
    n_r2_over_0.8_TX,
    delta_highLD_TX_minus_TT
  )
]

res <- merge(x, top_af, by = "Region_ID", all.x = TRUE, suffixes = c("", "_topAF"))
res <- merge(res, ld_sub, by = "Region_ID", all.x = TRUE)

final <- res[
  ,
  .(
    Region_ID,
    CHROM,
    START,
    END,
    Chromosome_type,
    Top_FST_SNP = Top_SNP,
    Max_SNP_FST,
    Top_deltaAF_SNP = SNP,
    AF_TT,
    AF_TX,
    delta_AF_TX_minus_TT,
    abs_delta_AF,
    delta_pi,
    TajimaD_TT,
    TajimaD_TX,
    delta_TajimaD,
    Ho_TT,
    Ho_TX,
    delta_Ho_TX_minus_TT,
    mean_r2_TT,
    mean_r2_TX,
    delta_mean_r2_TX_minus_TT,
    Nearest_gene = Gene_ID,
    Nearest_symbol = Gene_symbol,
    Gene_description = Description,
    Gene_biotype,
    Distance_to_gene_bp,
    Protein_coding_genes,
    Protein_coding_descriptions,
    Evidence_score_with_Het
  )
]

setorder(final, -Evidence_score_with_Het, -abs_delta_AF)

fwrite(final, out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Done\n")
cat("Output:", out_file, "\n")
print(final)
