#!/bin/bash
#SBATCH --job-name=LD50kb
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-67

set -euo pipefail

module load plink/1.90b6.24

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
BFILE=${BASE}/PreGWAS/result/GWAS_qc_named

LD=${BASE}/Genetics_Analysis/LD
WIN=${LD}/regions/genomewide_50kb_windows_for_LD.tsv
TMP=${LD}/plink/tmp_LD_windows
OUT=${LD}/results

mkdir -p ${TMP} ${OUT}

BLOCK=1000
START_LINE=$(( (${SLURM_ARRAY_TASK_ID} - 1) * ${BLOCK} + 2 ))
END_LINE=$(( ${SLURM_ARRAY_TASK_ID} * ${BLOCK} + 1 ))

OUTFILE=${OUT}/LD_50kb_summary_chunk_${SLURM_ARRAY_TASK_ID}.tsv

echo -e "CHROM\tBIN_START\tBIN_END\tN_SNPs\tGROUP\tN_PAIRS\tMEAN_R2\tN_R2_OVER_0.8" > ${OUTFILE}

sed -n "${START_LINE},${END_LINE}p" ${WIN} | while read CHROM BIN_START BIN_END N_SNPs
do
  for GROUP in TT TX
  do
    KEEP=${BASE}/Genetics_Analysis/Diversity/${GROUP}.keep
    PREFIX=${TMP}/${GROUP}_${CHROM}_${BIN_START}_${BIN_END}_${SLURM_ARRAY_TASK_ID}

    plink \
      --bfile ${BFILE} \
      --allow-extra-chr \
      --allow-no-sex \
      --keep ${KEEP} \
      --chr ${CHROM} \
      --from-bp ${BIN_START} \
      --to-bp ${BIN_END} \
      --r2 \
      --ld-window 999999 \
      --ld-window-kb 50 \
      --ld-window-r2 0 \
      --out ${PREFIX} >/dev/null 2>&1 || true

    if [ -s ${PREFIX}.ld ]; then
      awk -v chr=${CHROM} -v s=${BIN_START} -v e=${BIN_END} -v ns=${N_SNPs} -v g=${GROUP} '
        NR==1{
          for(i=1;i<=NF;i++) if($i=="R2") r2=i
          next
        }
        {
          n++
          sum+=$r2
          if($r2>=0.8) high++
        }
        END{
          if(n>0){print chr,s,e,ns,g,n,sum/n,high+0}
          else{print chr,s,e,ns,g,0,"NA",0}
        }
      ' OFS="\t" ${PREFIX}.ld >> ${OUTFILE}
    else
      echo -e "${CHROM}\t${BIN_START}\t${BIN_END}\t${N_SNPs}\t${GROUP}\t0\tNA\t0" >> ${OUTFILE}
    fi

    rm -f ${PREFIX}.ld ${PREFIX}.log ${PREFIX}.nosex
  done
done

echo "Done chunk ${SLURM_ARRAY_TASK_ID}"
