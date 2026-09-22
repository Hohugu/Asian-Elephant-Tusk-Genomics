#!/bin/bash
#SBATCH --job-name=camp_metrics
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/logs/%x_%j.err

set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction
ENV=$CAMP/envs/raisd_build

VCF=$CAMP/input/GWAS_qc_named_pheno.vcf.gz

export PATH=$ENV/bin:$CAMP/bin:/users/huguetho/vcftools/bin:$PATH
export LD_LIBRARY_PATH=$ENV/lib:${LD_LIBRARY_PATH:-}

mkdir -p $CAMP/metrics/fst
mkdir -p $CAMP/metrics/raisd
mkdir -p $CAMP/metrics/dxy_het

echo "Checking tools"
which bcftools
which tabix
which RAiSD
which vcftools

echo "Checking VCF"
ls -lh $VCF $VCF.tbi

bcftools query -l $VCF > $CAMP/input/GWAS_qc_named_pheno.samples.txt

echo "Sample count:"
wc -l $CAMP/input/GWAS_qc_named_pheno.samples.txt
wc -l $CAMP/input/TT.samples.vcf $CAMP/input/TX.samples.vcf

echo "1. FST TT vs TX"
cd $CAMP/metrics/fst

vcftools \
  --gzvcf $VCF \
  --weir-fst-pop $CAMP/input/TT.samples.vcf \
  --weir-fst-pop $CAMP/input/TX.samples.vcf \
  --out TT_vs_TX

echo "2. RAiSD muLD in TX"
cd $CAMP/metrics/raisd

RAiSD \
  -n TX_full_R \
  -I $VCF \
  -S $CAMP/input/TX.samples.vcf \
  -f \
  -R

awk 'BEGIN{OFS="\t"; print "CHROM","POS","window_start","window_end","muVAR","muSFS","muLD","mu"}
     /^\/\//{chrom=$2; next}
     NF>=7{print chrom,$1,$2,$3,$4,$5,$6,$7}' \
  RAiSD_Report.TX_full_R \
  > RAiSD_TX_muLD_full.tsv

echo "3. DXY and heterozygosity deviation"
cd $CAMP/metrics/dxy_het

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
  $VCF | \
python3 $CAMP/scripts/compute_dxy_hetdev_stream.py \
  --sample-order $CAMP/input/GWAS_qc_named_pheno.samples.txt \
  --tt $CAMP/input/TT.samples.vcf \
  --tx $CAMP/input/TX.samples.vcf \
  --out DXY_HetDev_TX.full.tsv.gz

echo "Done"
ls -lh \
  $CAMP/metrics/fst/TT_vs_TX.weir.fst \
  $CAMP/metrics/raisd/RAiSD_Report.TX_full_R \
  $CAMP/metrics/raisd/RAiSD_TX_muLD_full.tsv \
  $CAMP/metrics/dxy_het/DXY_HetDev_TX.full.tsv.gz
