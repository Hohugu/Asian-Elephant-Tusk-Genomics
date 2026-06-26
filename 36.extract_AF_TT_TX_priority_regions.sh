#!/bin/bash

set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
FA=${BASE}/Genetics_Analysis/Functional_annotation
VCFDIR=${FA}/vcf
OUT=${FA}/tables

mkdir -p ${OUT}

for REGION in FST_region_18 FST_region_9 FST_region_56
do
  echo "Processing ${REGION}"

  bcftools view \
    -S ${BASE}/Genetics_Analysis/Diversity/TT.vcf.samples \
    ${VCFDIR}/${REGION}.vcf.gz \
  | bcftools +fill-tags \
    -- -t AF,AC,AN \
  | bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
    > ${OUT}/${REGION}_TT_AF.tsv

  bcftools view \
    -S ${BASE}/Genetics_Analysis/Diversity/TX.vcf.samples \
    ${VCFDIR}/${REGION}.vcf.gz \
  | bcftools +fill-tags \
    -- -t AF,AC,AN \
  | bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
    > ${OUT}/${REGION}_TX_AF.tsv

done

echo "Done"
