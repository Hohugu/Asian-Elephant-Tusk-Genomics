#!/bin/bash
#SBATCH --job-name=campbell_pc3
#SBATCH --output=campbell_pc3.out
#SBATCH --error=campbell_pc3.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

set -euo pipefail

module load bedtools

# =========================
# Input files
# =========================

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"

GFF="/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

GWAS="${BASE}/GWAS/result/GEMMA_PC3/output/tusk_GEMMA_LMM_PC3.assoc.txt"

MAP="${BASE}/GWAS/tables_GEMMA_male/CM_to_NC.map"

OUTDIR="${BASE}/GWAS/tables_GEMMA_PC3/Campbell_candidates_final"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# =========================
# 1. Extract Campbell candidate genes from GFF
# =========================
# Candidate genes:
# AMELX / amelogenin
# MEP1A / meprin A subunit alpha
# ENAM / enamelin
# AMBN / ameloblastin
# AMTN / amelotin
# ODAM / odontogenic ameloblast-associated

awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene" &&
$9 ~ /amelogenin|meprin A subunit alpha|enamelin|ameloblastin|amelotin|odontogenic/ {

    gene="NA"

    if($9 ~ /meprin A subunit alpha/) gene="MEP1A"
    else if($9 ~ /amelogenin X-linked/) gene="AMELX"
    else if($9 ~ /amelogenin%2C X isoform/) gene="AMELX_like"
    else if($9 ~ /enamelin/) gene="ENAM"
    else if($9 ~ /ameloblastin/) gene="AMBN"
    else if($9 ~ /amelotin/) gene="AMTN"
    else if($9 ~ /odontogenic/) gene="ODAM"

    loc="NA"
    if(match($9,/gene=([^;]+)/,a)) loc=a[1]

    desc="NA"
    if(match($9,/description=([^;]+)/,b)) desc=b[1]

    print $1,$4-1,$5,gene"_"loc,desc
}
' "$GFF" | sort -k1,1 -k2,2n \
> Campbell_candidate_genes.bed

# =========================
# 2. Create 500 kb, 1 Mb and 2 Mb windows
# =========================

for W in 500000 1000000 2000000
do
    awk -v W="$W" 'BEGIN{OFS="\t"}{
        s=$2-W
        if(s<0)s=0
        e=$3+W
        print $1,s,e,$4,$5
    }' Campbell_candidate_genes.bed \
    > Campbell_candidate_genes_${W}.bed
done

# =========================
# 3. Convert GEMMA results from CM IDs to NC IDs
# =========================
# Important:
# CM044048.1 corresponds to NC_064847.1, which was absent from the old map.
# We add it manually as a fallback.

awk -v MAP="$MAP" '
BEGIN{
    FS=OFS="\t"
    while((getline<MAP)>0){
        map[$1]=$2
    }

    # manual correction / fallback
    map["CM044048.1"]="NC_064847.1"
}
NR==1{next}
{
    chr=$1
    if(chr in map) chr=map[chr]

    # GEMMA columns:
    # chr=$1, rs=$2, ps=$3, p_wald=$13
    print chr,$3-1,$3,$2,$13
}
' "$GWAS" \
> GEMMA_PC3_allSNPs.NC.bed

# =========================
# 4. Compute Bonferroni threshold
# =========================

NTESTS=$(awk 'NR>1{n++} END{print n}' "$GWAS")
BONF=$(awk -v n="$NTESTS" 'BEGIN{print 0.05/n}')

echo -e "N_tests\tBonferroni_threshold" > GEMMA_PC3_Bonferroni_threshold.tsv
echo -e "${NTESTS}\t${BONF}" >> GEMMA_PC3_Bonferroni_threshold.tsv

# =========================
# 5. Extract Bonferroni-significant SNPs
# =========================

awk -v MAP="$MAP" -v BONF="$BONF" '
BEGIN{
    FS=OFS="\t"
    while((getline<MAP)>0){
        map[$1]=$2
    }

    # manual correction / fallback
    map["CM044048.1"]="NC_064847.1"
}
NR==1{next}
$13 < BONF {
    chr=$1
    if(chr in map) chr=map[chr]

    print chr,$3-1,$3,$2,$13
}
' "$GWAS" \
> GEMMA_PC3_Bonferroni_hits.NC.bed

# =========================
# 6. Retrieve all SNP p-values inside and around candidate genes
# =========================

bedtools intersect \
-a GEMMA_PC3_allSNPs.NC.bed \
-b Campbell_candidate_genes.bed \
-wa -wb \
> GEMMA_PC3_SNPs_inside_Campbell_candidates.tsv

for W in 500000 1000000 2000000
do
    bedtools intersect \
    -a GEMMA_PC3_allSNPs.NC.bed \
    -b Campbell_candidate_genes_${W}.bed \
    -wa -wb \
    > GEMMA_PC3_SNPs_${W}_Campbell_candidates.tsv
done

# =========================
# 7. Best p-value per candidate gene / region
# =========================

for REGION in inside 500000 1000000 2000000
do
    if [ "$REGION" = "inside" ]; then
        INFILE="GEMMA_PC3_SNPs_inside_Campbell_candidates.tsv"
    else
        INFILE="GEMMA_PC3_SNPs_${REGION}_Campbell_candidates.tsv"
    fi

    awk -F'\t' '
    BEGIN{
        OFS="\t"
        print "Candidate_gene","Best_SNP","Best_p_wald"
    }
    {
        gene=$9
        snp=$4
        p=$5

        if(!(gene in best) || p < best[gene]){
            best[gene]=p
            bestsnp[gene]=snp
        }
    }
    END{
        for(g in best){
            print g,bestsnp[g],best[g]
        }
    }
    ' "$INFILE" | sort -k3,3g \
    > GEMMA_PC3_best_pvalues_${REGION}_Campbell_candidates.tsv
done

# =========================
# 8. Bonferroni-significant SNPs near candidate genes
# =========================

bedtools intersect \
-a GEMMA_PC3_Bonferroni_hits.NC.bed \
-b Campbell_candidate_genes.bed \
-wa -wb \
> GEMMA_PC3_Bonferroni_hits_inside_Campbell_candidates.tsv

for W in 500000 1000000 2000000
do
    bedtools intersect \
    -a GEMMA_PC3_Bonferroni_hits.NC.bed \
    -b Campbell_candidate_genes_${W}.bed \
    -wa -wb \
    > GEMMA_PC3_Bonferroni_hits_${W}_Campbell_candidates.tsv
done

# =========================
# 9. Final summary
# =========================

{
echo "Candidate gene positions:"
cat Campbell_candidate_genes.bed
echo
echo "Bonferroni threshold:"
cat GEMMA_PC3_Bonferroni_threshold.tsv
echo
echo "Best p-values inside genes:"
cat GEMMA_PC3_best_pvalues_inside_Campbell_candidates.tsv
echo
echo "Best p-values within 500 kb:"
cat GEMMA_PC3_best_pvalues_500000_Campbell_candidates.tsv
echo
echo "Best p-values within 1 Mb:"
cat GEMMA_PC3_best_pvalues_1000000_Campbell_candidates.tsv
echo
echo "Best p-values within 2 Mb:"
cat GEMMA_PC3_best_pvalues_2000000_Campbell_candidates.tsv
echo
echo "Bonferroni hits inside candidate genes:"
wc -l GEMMA_PC3_Bonferroni_hits_inside_Campbell_candidates.tsv
echo
echo "Bonferroni hits within 500 kb:"
wc -l GEMMA_PC3_Bonferroni_hits_500000_Campbell_candidates.tsv
echo
echo "Bonferroni hits within 1 Mb:"
wc -l GEMMA_PC3_Bonferroni_hits_1000000_Campbell_candidates.tsv
echo
echo "Bonferroni hits within 2 Mb:"
wc -l GEMMA_PC3_Bonferroni_hits_2000000_Campbell_candidates.tsv
} > Campbell_candidates_GEMMA_PC3_summary.txt

echo "Done."
echo "Output directory: $OUTDIR"
