#!/bin/bash
#SBATCH --job-name=annot_bonf
#SBATCH --output=annot_bonf.out
#SBATCH --error=annot_bonf.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

set -euo pipefail

module load bedtools

BASE="/scratch/project_2000886/Hoedric/GWAS_2025"
GFF="/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

OUTDIR="${BASE}/GWAS/tables_GEMMA_PC3/Campbell_candidates_final"
cd "$OUTDIR"

# =========================
# 1. Create BED file of all genes from reference GFF
# =========================

awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene" {
    gene="NA"
    desc="NA"
    biotype="NA"

    if(match($9,/Name=([^;]+)/,a)) gene=a[1]
    if(match($9,/description=([^;]+)/,b)) desc=b[1]
    if(match($9,/gene_biotype=([^;]+)/,c)) biotype=c[1]

    print $1,$4-1,$5,gene,desc,biotype
}
' "$GFF" | sort -k1,1 -k2,2n \
> all_genes_from_GFF.sorted.bed

# =========================
# 2. Annotate Bonferroni SNPs within 500 kb, 1 Mb and 2 Mb
# =========================

for W in 500000 1000000 2000000
do
    INFILE="GEMMA_PC3_Bonferroni_hits_${W}_Campbell_candidates.tsv"

    if [ ! -s "$INFILE" ]; then
        echo -e "SNP_chr\tSNP_pos\tSNP\tp_wald\tCandidate_gene\tDistance_to_candidate_gene\tOverlapping_gene\tOverlapping_gene_description\tOverlapping_gene_type" \
        > GEMMA_PC3_Bonferroni_hits_${W}_annotated.tsv
        continue
    fi

    # Extract SNPs only
    awk -F'\t' 'BEGIN{OFS="\t"}{
        print $1,$2,$3,$4,$5
    }' "$INFILE" \
    | sort -k1,1 -k2,2n \
    > tmp_hits_${W}.bed

    # Find overlapping genes
    bedtools intersect \
    -a tmp_hits_${W}.bed \
    -b all_genes_from_GFF.sorted.bed \
    -wa -wb \
    > tmp_hits_${W}_overlap.tsv || true

    # Find nearest genes if no overlap
    bedtools closest \
    -a tmp_hits_${W}.bed \
    -b all_genes_from_GFF.sorted.bed \
    -d \
    > tmp_hits_${W}_nearest.tsv

    # Build final table with candidate gene and distance to candidate
    awk -F'\t' '
    BEGIN{
        OFS="\t"
        print "SNP_chr","SNP_pos","SNP","p_wald","Candidate_gene","Distance_to_candidate_gene","Overlapping_or_nearest_gene","Gene_description","Gene_type","Distance_to_gene"
    }

    FNR==NR{
        # overlap file
        key=$1":"$3":"$4
        ogene[key]=$9
        odesc[key]=$10
        otype[key]=$11
        odist[key]=0
        next
    }

    {
        # nearest file
        key=$1":"$3":"$4

        snp_chr=$1
        snp_pos=$3
        snp=$4
        p=$5

        gene_chr=$6
        gene_start=$7
        gene_end=$8
        ngene=$9
        ndesc=$10
        ntype=$11
        ndist=$12

        if(key in ogene){
            final_gene=ogene[key]
            final_desc=odesc[key]
            final_type=otype[key]
            final_dist=0
        } else {
            final_gene=ngene
            final_desc=ndesc
            final_type=ntype
            final_dist=ndist
        }

        print snp_chr,snp_pos,snp,p,"NA","NA",final_gene,final_desc,final_type,final_dist
    }
    ' tmp_hits_${W}_overlap.tsv tmp_hits_${W}_nearest.tsv \
    > GEMMA_PC3_Bonferroni_hits_${W}_annotated.tsv

done

# =========================
# 3. Specific AMELX lead SNP summary
# =========================
# This captures the already observed Bonferroni SNP near AMELX.

grep "AMELX_LOC126069472" GEMMA_PC3_Bonferroni_hits_2000000_Campbell_candidates.tsv \
> AMELX_2Mb_Bonferroni_hits.tsv || true

if [ -s AMELX_2Mb_Bonferroni_hits.tsv ]; then

awk -F'\t' '
BEGIN{
    OFS="\t"
    print "SNP_chr","SNP_pos","SNP","p_wald","AMELX_start","AMELX_end","Distance_to_AMELX"
}
{
    snp_pos=$3
    amelx_start=168725165
    amelx_end=168729418

    if(snp_pos < amelx_start) dist=amelx_start-snp_pos
    else if(snp_pos > amelx_end) dist=snp_pos-amelx_end
    else dist=0

    print $1,snp_pos,$4,$5,amelx_start,amelx_end,dist
}
' AMELX_2Mb_Bonferroni_hits.tsv \
> AMELX_2Mb_Bonferroni_hits_distance.tsv

fi

# =========================
# 4. Genes between AMELX and its lead SNP
# =========================

if [ -s AMELX_2Mb_Bonferroni_hits_distance.tsv ]; then

LEAD_POS=$(awk 'NR==2{print $2}' AMELX_2Mb_Bonferroni_hits_distance.tsv)

awk -F'\t' -v LEAD="$LEAD_POS" '
BEGIN{OFS="\t"}
$1=="NC_064846.1" &&
$3=="gene" &&
$5>=168729418 &&
$4<=LEAD {
    gene="NA"
    desc="NA"
    biotype="NA"

    if(match($9,/Name=([^;]+)/,a)) gene=a[1]
    if(match($9,/description=([^;]+)/,b)) desc=b[1]
    if(match($9,/gene_biotype=([^;]+)/,c)) biotype=c[1]

    print $1,$4,$5,gene,desc,biotype
}
' "$GFF" \
> genes_between_AMELX_and_leadSNP.tsv

fi

# =========================
# 5. Final report
# =========================

{
echo "Annotated Bonferroni hits within 500 kb:"
cat GEMMA_PC3_Bonferroni_hits_500000_annotated.tsv
echo
echo "Annotated Bonferroni hits within 1 Mb:"
cat GEMMA_PC3_Bonferroni_hits_1000000_annotated.tsv
echo
echo "Annotated Bonferroni hits within 2 Mb:"
cat GEMMA_PC3_Bonferroni_hits_2000000_annotated.tsv
echo
echo "AMELX 2 Mb lead SNP distance:"
if [ -s AMELX_2Mb_Bonferroni_hits_distance.tsv ]; then
    cat AMELX_2Mb_Bonferroni_hits_distance.tsv
else
    echo "No AMELX 2 Mb Bonferroni hit found"
fi
echo
echo "Genes between AMELX and lead SNP:"
if [ -s genes_between_AMELX_and_leadSNP.tsv ]; then
    cat genes_between_AMELX_and_leadSNP.tsv
else
    echo "No genes listed"
fi
} > GEMMA_PC3_Campbell_Bonferroni_annotation_summary.txt

echo "Done."
echo "Main summary:"
echo "${OUTDIR}/GEMMA_PC3_Campbell_Bonferroni_annotation_summary.txt"
