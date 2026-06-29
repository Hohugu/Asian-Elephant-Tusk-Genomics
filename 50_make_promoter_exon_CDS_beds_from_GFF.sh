#!/bin/bash

set -euo pipefail

GFF="/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

OUTDIR="/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Functional_annotation/regulatory"
mkdir -p "${OUTDIR}"

# ------------------------------------------------------------------
# Genes
# ------------------------------------------------------------------

awk '
BEGIN{OFS="\t"}
$1 !~ /^#/ && $3=="gene" {
  print $1,$4-1,$5
}
' "$GFF" \
| sort -k1,1 -k2,2n \
> "${OUTDIR}/genes.bed"

# ------------------------------------------------------------------
# Exons
# ------------------------------------------------------------------

awk '
BEGIN{OFS="\t"}
$1 !~ /^#/ && $3=="exon" {
  print $1,$4-1,$5
}
' "$GFF" \
| sort -k1,1 -k2,2n \
> "${OUTDIR}/exons.bed"

# ------------------------------------------------------------------
# CDS
# ------------------------------------------------------------------

awk '
BEGIN{OFS="\t"}
$1 !~ /^#/ && $3=="CDS" {
  print $1,$4-1,$5
}
' "$GFF" \
| sort -k1,1 -k2,2n \
> "${OUTDIR}/CDS.bed"

# ------------------------------------------------------------------
# Promoters: -2000/+500 autour du TSS
# ------------------------------------------------------------------

awk '
BEGIN{OFS="\t"}

$1 !~ /^#/ && $3=="gene" {

  strand=$7

  if(strand=="+"){
      start=$4-2000
      end=$4+500
  } else {
      start=$5-500
      end=$5+2000
  }

  if(start<0) start=0

  print $1,start,end
}
' "$GFF" \
| sort -k1,1 -k2,2n \
> "${OUTDIR}/promoters_2kb.bed"

echo "Done"
echo "${OUTDIR}"
wc -l "${OUTDIR}"/*.bed
