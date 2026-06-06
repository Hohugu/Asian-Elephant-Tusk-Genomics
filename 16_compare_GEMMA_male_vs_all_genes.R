#!/usr/bin/env Rscript

library(data.table)

base <- "/scratch/project_2000886/Hoedric/GWAS_2025/GWAS"

gff_file <- "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

all_hits_file  <- file.path(base, "result/GEMMA/output/GEMMA_Bonferroni_hits.txt")
male_hits_file <- file.path(base, "tables_GEMMA_male/GEMMA_male_Bonferroni_hits.tsv")

outdir <- file.path(base, "tables_GEMMA_compare")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Checking existing GEMMA annotation files...\n")
print(list.files(file.path(base, "result/GEMMA/output"), pattern="annotation|annot|gene|hits", full.names=TRUE))
print(list.files(file.path(base, "tables_GEMMA"), pattern="annotation|annot|gene|hits", full.names=TRUE))

# -----------------------------
# CM -> NC mapping
# -----------------------------
cm <- sprintf("CM044%03d.1", 20:47)
nc <- sprintf("NC_064%d.1", 819:846)
chr_map <- data.table(chr=cm, chr_nc=nc)

# -----------------------------
# Read GEMMA hits
# -----------------------------
all_hits <- fread(all_hits_file, header=FALSE)
setnames(all_hits, c(
  "chr","rs","ps","n_miss","allele1","allele0","af","beta","se",
  "logl_H1","l_remle","l_mle","p_wald","p_lrt","p_score"
))

male_hits <- fread(male_hits_file)

all_hits[, analysis := "GEMMA_all"]
male_hits[, analysis := "GEMMA_male"]

all_hits <- merge(all_hits, chr_map, by="chr", all.x=TRUE)
male_hits <- merge(male_hits, chr_map, by="chr", all.x=TRUE)

# -----------------------------
# Parse GFF genes
# -----------------------------
gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(gff_file)),
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE
)

gff <- gff[V3 == "gene"]
setnames(gff, c("seqid","source","type","start","end","score","strand","phase","attributes"))

extract_attr <- function(x, key) {
  m <- regmatches(x, regexpr(paste0(key, "=[^;]+"), x))
  sub(paste0(key, "="), "", m)
}

gff[, gene_id := extract_attr(attributes, "ID")]
gff[, gene_name := extract_attr(attributes, "gene")]
gff[is.na(gene_name) | gene_name == "", gene_name := extract_attr(attributes, "Name")]
gff[is.na(gene_name) | gene_name == "", gene_name := gene_id]

genes <- gff[, .(
  gene_chr = seqid,
  gene_start = start - 1,
  gene_end = end,
  gene_id,
  gene_name,
  strand
)]

# -----------------------------
# Nearest gene function
# -----------------------------
nearest_gene <- function(hits) {
  res <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
    h <- hits[i]
    g <- genes[gene_chr == h$chr_nc]
    
    if (nrow(g) == 0) {
      return(data.table(
        analysis=h$analysis, chr=h$chr, chr_nc=h$chr_nc, rs=h$rs,
        snp_pos=h$ps, p_wald=h$p_wald,
        gene_chr=NA, gene_start=NA, gene_end=NA,
        gene_id=NA, gene_name=NA, strand=NA, distance_bp=NA
      ))
    }
    
    g[, distance_bp := fifelse(
      h$ps >= gene_start & h$ps <= gene_end,
      0,
      pmin(abs(h$ps - gene_start), abs(h$ps - gene_end))
    )]
    
    g <- g[order(distance_bp)][1]
    
    data.table(
      analysis=h$analysis,
      chr=h$chr,
      chr_nc=h$chr_nc,
      rs=h$rs,
      snp_pos=h$ps,
      p_wald=h$p_wald,
      gene_chr=g$gene_chr,
      gene_start=g$gene_start,
      gene_end=g$gene_end,
      gene_id=g$gene_id,
      gene_name=g$gene_name,
      strand=g$strand,
      distance_bp=g$distance_bp
    )
  }))
  
  return(res)
}

cat("Annotating all GEMMA hits...\n")
all_annot <- nearest_gene(all_hits)

cat("Annotating male-only GEMMA hits...\n")
male_annot <- nearest_gene(male_hits)

# -----------------------------
# Write per-analysis tables
# -----------------------------
fwrite(
  all_annot,
  file.path(outdir, "GEMMA_all_Bonferroni_hits_nearest_genes.tsv"),
  sep="\t"
)

fwrite(
  male_annot,
  file.path(outdir, "GEMMA_male_Bonferroni_hits_nearest_genes.tsv"),
  sep="\t"
)

# -----------------------------
# Shared genes
# -----------------------------
all_genes <- unique(all_annot[!is.na(gene_name), .(gene_name)])
male_genes <- unique(male_annot[!is.na(gene_name), .(gene_name)])

shared_genes <- merge(all_genes, male_genes, by="gene_name")

shared_table <- merge(
  all_annot,
  shared_genes,
  by="gene_name"
)

shared_table <- merge(
  shared_table,
  male_annot,
  by="gene_name",
  suffixes=c("_all","_male"),
  allow.cartesian=TRUE
)

fwrite(
  shared_table,
  file.path(outdir, "GEMMA_shared_genes_all_vs_male.tsv"),
  sep="\t"
)

# -----------------------------
# SNP overlap
# -----------------------------
shared_snps <- intersect(all_annot$rs, male_annot$rs)

fwrite(
  data.table(shared_snps=shared_snps),
  file.path(outdir, "GEMMA_shared_SNPs_all_vs_male.tsv"),
  sep="\t"
)

cat("Done.\n")
cat("All hits:", nrow(all_annot), "\n")
cat("Male hits:", nrow(male_annot), "\n")
cat("Shared SNPs:", length(shared_snps), "\n")
cat("Shared genes:", nrow(shared_genes), "\n")
cat("Output dir:", outdir, "\n")
