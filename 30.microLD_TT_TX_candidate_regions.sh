#!/bin/bash

set -euo pipefail

module load plink/1.90b6.24

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
BFILE=${BASE}/PreGWAS/result/GWAS_qc_named

LD=${BASE}/Genetics_Analysis/LD
REG=${LD}/regions
OUT=${LD}/results
PLINKOUT=${LD}/plink

mkdir -p ${REG} ${OUT} ${PLINKOUT}

# Candidate regions for local LD / microLD
cat > ${REG}/candidate_regions_microLD.tsv <<EOF
Region	CHR	Start	End	Reason
FST_region_18	CM044020.1	77550001	78100000	FST_pi_TajimaD_ZNF184_ZNF391
FST_region_9	CM044022.1	120450001	121000000	FST_pi_TajimaD_U6
FST_region_56	CM044021.1	136800001	137350000	FST_pi_TajimaD_intergenic
PDK3	CM044047.1	154200000	154700000	GWAS_mixed_FST
AMELX	CM044047.1	168300000	170200000	Campbell_tooth_gene
MEP1A_PLA2G7	CM044020.1	110000000	112000000	Campbell_tooth_region
MEP1B	CM044030.1	38500000	41500000	Campbell_tooth_gene
EOF

echo "Running local LD for TT and TX"

tail -n +2 ${REG}/candidate_regions_microLD.tsv | while read REGION CHR START END REASON
do
  echo "======================================"
  echo ${REGION} ${CHR}:${START}-${END}
  echo "======================================"

  plink \
    --bfile ${BFILE} \
    --allow-extra-chr \
    --allow-no-sex \
    --keep ${BASE}/Genetics_Analysis/Diversity/TT.keep \
    --chr ${CHR} \
    --from-bp ${START} \
    --to-bp ${END} \
    --r2 \
    --ld-window 999999 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --out ${PLINKOUT}/${REGION}_TT_LD

  plink \
    --bfile ${BFILE} \
    --allow-extra-chr \
    --allow-no-sex \
    --keep ${BASE}/Genetics_Analysis/Diversity/TX.keep \
    --chr ${CHR} \
    --from-bp ${START} \
    --to-bp ${END} \
    --r2 \
    --ld-window 999999 \
    --ld-window-kb 1000 \
    --ld-window-r2 0 \
    --out ${PLINKOUT}/${REGION}_TX_LD

done

echo "Done"
echo "LD outputs in: ${PLINKOUT}"
