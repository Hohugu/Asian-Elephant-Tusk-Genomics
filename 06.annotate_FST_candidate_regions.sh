#!/bin/bash

set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Selection/FST

REGIONS=${BASE}/tables/TT_vs_TX_FST_50kb_candidate_regions_merged.tsv

GFF=/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff
CM_NC=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/tables_GEMMA_male/CM_to_NC.map

TMP_CM=${BASE}/tables/fst_regions.CM.bed
TMP_NC=${BASE}/tables/fst_regions.NC.bed
GENES=${BASE}/tables/all_genes.bed

awk 'BEGIN{FS=OFS="\t"}
NR>1{
  print $2,$3,$4,$1
}' "$REGIONS" > "$TMP_CM"

echo "Converting FST regions CM -> NC"

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{
  map[$1]=$2
  next
}
{
  if($1 in map){
    $1=map[$1]
    print
  }
}' "$CM_NC" "$TMP_CM" > "$TMP_NC"

echo "Extracting genes from GFF"

awk '
BEGIN{FS=OFS="\t"}
$3=="gene"{

  chr=$1
  start=$4-1
  end=$5

  id="NA"
  name="NA"
  desc="NA"

  if(match($9,/ID=gene-([^;]+)/)){
    id=substr($9,RSTART+8,RLENGTH-8)
  }

  if(match($9,/gene=([^;]+)/)){
    name=substr($9,RSTART+5,RLENGTH-5)
  } else {
    name=id
  }

  if(match($9,/description=([^;]+)/)){
    desc=substr($9,RSTART+12,RLENGTH-12)
  }

  print chr,start,end,id,name,desc
}
' "$GFF" > "$GENES"

echo "Intersect"

bedtools intersect \
-wa -wb \
-a "$TMP_NC" \
-b "$GENES" \
> ${BASE}/tables/FST_regions_gene_overlap.tsv

echo "Done"
