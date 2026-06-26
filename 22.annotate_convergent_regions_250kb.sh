#!/bin/bash

set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
DIV=${BASE}/Genetics_Analysis/Diversity
OUT=${DIV}/convergent_regions_annotation_250kb

mkdir -p ${OUT}

GFF=/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff
MAP=${BASE}/GWAS/tables_GEMMA_male/CM_to_NC.map
IN=${DIV}/FST_pi_TajimaD_convergence.tsv

# regions ±250 kb, CM coordinates
awk 'BEGIN{FS=OFS="\t"}
NR>1{
  s=$4-250000
  if(s<1){s=1}
  e=$5+250000
  print $3,s,e,$2,$1,$4,$5
}' ${IN} | sort -u > ${OUT}/convergent_regions_250kb.CM.bed

# convert CM -> NC
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{map[$1]=$2; next}
($1 in map){
  $1=map[$1]
  print
}' ${MAP} ${OUT}/convergent_regions_250kb.CM.bed \
> ${OUT}/convergent_regions_250kb.NC.bed

# extract genes from GFF
awk 'BEGIN{FS=OFS="\t"}
$3=="gene"{
  id="NA"; symbol="NA"; desc="NA"; biotype="NA"

  if(match($9,/ID=gene-([^;]+)/)){
    id=substr($9,RSTART+8,RLENGTH-8)
  }
  if(match($9,/gene=([^;]+)/)){
    symbol=substr($9,RSTART+5,RLENGTH-5)
  } else {
    symbol=id
  }
  if(match($9,/description=([^;]+)/)){
    desc=substr($9,RSTART+12,RLENGTH-12)
  }
  if(match($9,/gene_biotype=([^;]+)/)){
    biotype=substr($9,RSTART+13,RLENGTH-13)
  }

  print $1,$4-1,$5,id,symbol,desc,biotype
}' ${GFF} > ${OUT}/all_genes.bed

bedtools intersect \
  -wa -wb \
  -a ${OUT}/convergent_regions_250kb.NC.bed \
  -b ${OUT}/all_genes.bed \
  > ${OUT}/convergent_regions_250kb_gene_overlap.tsv

echo "Done"
echo "Output:"
echo "${OUT}/convergent_regions_250kb_gene_overlap.tsv"
wc -l ${OUT}/convergent_regions_250kb_gene_overlap.tsv
