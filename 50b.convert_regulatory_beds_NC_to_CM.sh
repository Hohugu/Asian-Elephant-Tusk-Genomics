#!/bin/bash
set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
REG=$BASE/Genetics_Analysis/Functional_annotation/regulatory
MAP=$BASE/GWAS/tables_GEMMA_male/CM_to_NC.map

for f in genes exons CDS promoters_2kb
do
  awk 'BEGIN{OFS="\t"} NR==FNR{cm[$2]=$1; next} ($1 in cm){$1=cm[$1]; print}' \
    "$MAP" "$REG/${f}.bed" \
  | sort -k1,1 -k2,2n \
  > "$REG/${f}.CM.bed"
done

echo "Done"
wc -l "$REG"/*.CM.bed
(base) [huguetho@puhti-login11 scripts]$ head 52b_convert_regulatory_beds_NC_to_CM.sh -n 10000
#!/bin/bash
set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
REG=$BASE/Genetics_Analysis/Functional_annotation/regulatory
MAP=$BASE/GWAS/tables_GEMMA_male/CM_to_NC.map

for f in genes exons CDS promoters_2kb
do
  awk 'BEGIN{OFS="\t"} NR==FNR{cm[$2]=$1; next} ($1 in cm){$1=cm[$1]; print}' \
    "$MAP" "$REG/${f}.bed" \
  | sort -k1,1 -k2,2n \
  > "$REG/${f}.CM.bed"
done

echo "Done"
wc -l "$REG"/*.CM.bed
