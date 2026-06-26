#!/bin/bash

set -euo pipefail

module load biokit

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
VCF=${BASE}/Genetics_Analysis/Diversity/GWAS_qc_named.vcf.gz

OUT=${BASE}/Genetics_Analysis/Functional_annotation
mkdir -p ${OUT}/{regions,vcf,tables}

cat > ${OUT}/regions/priority_regions.bed <<EOF
CM044020.1	77800000	77850000	FST_region_18
CM044022.1	120700000	120750000	FST_region_9
CM044021.1	137050000	137100000	FST_region_56
EOF

while read CHR START END REGION
do
  echo "Extracting ${REGION} ${CHR}:${START}-${END}"

  bcftools view \
    -r ${CHR}:$((START+1))-${END} \
    -Oz \
    -o ${OUT}/vcf/${REGION}.vcf.gz \
    ${VCF}

  tabix -f -p vcf ${OUT}/vcf/${REGION}.vcf.gz

  bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\n' \
    ${OUT}/vcf/${REGION}.vcf.gz \
    > ${OUT}/tables/${REGION}_variants.tsv

done < ${OUT}/regions/priority_regions.bed

echo "Done"
