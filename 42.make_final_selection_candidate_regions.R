suppressPackageStartupMessages({
  library(data.table)
})

base <- Sys.getenv("BASE", "/scratch/project_2000886/Hoedric/GWAS_2025")
final <- file.path(base, "Genetics_Analysis/Final_candidate_loci")
dir.create(file.path(final, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(final, "notes"), recursive = TRUE, showWarnings = FALSE)

src_final <- file.path(final, "tables", "integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv")
src_original <- file.path(base, "Genetics_Analysis/Diversity/integrated_candidate_regions_FST_pi_TajimaD_LD_Het.tsv")

if (!file.exists(src_final)) {
  if (!file.exists(src_original)) {
    stop("Missing source file: ", src_original)
  }
  file.copy(src_original, src_final, overwrite = TRUE)
}

regions <- fread(src_final)

keep <- c(
  "Region_ID", "CHROM", "START", "END", "Chromosome_type",
  "Mean_FST", "Max_SNP_FST", "Top_SNP", "Top_SNP_POS",
  "delta_pi", "TajimaD_TT", "TajimaD_TX", "delta_TajimaD",
  "Ho_TT", "Ho_TX", "delta_Ho_TX_minus_TT", "Het_direction",
  "Evidence_score_with_Het", "Protein_coding_genes", "Protein_coding_descriptions", "Direction"
)
keep <- keep[keep %in% names(regions)]
sel <- regions[, ..keep]

sel[, Final_category := "Selection/differentiation candidate region"]
sel[, Interpretation := "Differentiated region with reduced TX diversity/heterozygosity; compatible with selection or haplotypic differentiation"]
sel[, Caution := "Not proof of selection; drift, population structure and technical artifacts are not fully excluded"]

out <- file.path(final, "tables", "Final_selection_candidate_regions.tsv")
fwrite(sel, out, sep = "\t")

note <- file.path(final, "notes", "Final_selection_candidate_regions_summary.txt")
cat(
  "Final selection/differentiation candidate regions\n",
  "================================================\n\n",
  "These regions are separate from the Tier 1/2/3 GWAS candidate-locus table.\n",
  "They summarize genomic differentiation between TT and TX using FST, pi, TajimaD, heterozygosity and LD evidence.\n\n",
  "Number of regions: ", nrow(sel), "\n\n",
  "Interpretation: compatible with selection or haplotypic differentiation, but drift/population structure is not excluded.\n",
  file = note,
  sep = ""
)

cat("Done\n")
cat("Selection regions table:", out, "\n")
cat("Summary note:", note, "\n\n")
cat("Regions:\n")
print(sel[, .(Region_ID, CHROM, START, END, Top_SNP, Max_SNP_FST, delta_pi, delta_TajimaD, Het_direction, Evidence_score_with_Het)])
