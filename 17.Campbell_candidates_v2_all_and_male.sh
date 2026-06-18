#!/bin/bash
#SBATCH --job-name=campbell_v2
#SBATCH --output=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/campbell_v2.out
#SBATCH --error=/scratch/project_2000886/Hoedric/GWAS_2025/GWAS/logs/campbell_v2.err
#SBATCH --account=project_2000886
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail

BASE="/scratch/project_2000886/Hoedric/GWAS_2025/GWAS"

CAND="${BASE}/script/Campbell_candidates_v2_coordinates.tsv"

cat > "$CAND" << EOF
Gene	CHR	Start	End	Description
AMBN	CM044024.1	87354290	87366652	ameloblastin
AMTN	CM044024.1	87416785	87432258	amelotin
ENAM	CM044024.1	87316884	87332008	enamelin
ODAM	CM044024.1	87832910	87841335	odontogenic ameloblast associated
AMELX	CM044047.1	168725166	168729418	amelogenin X-linked
MEP1A	CM044020.1	111679948	111709525	meprin A subunit alpha
MEP1B	CM044030.1	39647956	39684612	meprin A subunit beta
PLA2G7	CM044020.1	111553602	111601754	phospholipase A2 group VII
EOF

run_analysis () {

    NAME="$1"
    GWAS="$2"
    OUTDIR="$3"

    mkdir -p "$OUTDIR"

    echo "======================================"
    echo "Running $NAME"
    echo "GWAS file: $GWAS"
    echo "Output dir: $OUTDIR"
    echo "======================================"

    NTESTS=$(awk 'NR>1{n++} END{print n}' "$GWAS")
    BONF=$(awk -v n="$NTESTS" 'BEGIN{print 0.05/n}')

    echo -e "Analysis\tN_tests\tBonferroni_threshold" > "${OUTDIR}/${NAME}_Bonferroni_threshold.tsv"
    echo -e "${NAME}\t${NTESTS}\t${BONF}" >> "${OUTDIR}/${NAME}_Bonferroni_threshold.tsv"

    awk -v CAND="$CAND" -v BONF="$BONF" -v OUT="$OUTDIR" -v NAME="$NAME" '
    BEGIN{
        FS=OFS="\t"

        while((getline < CAND) > 0){
            if($1=="Gene") continue
            n++
            gene[n]=$1
            chr[n]=$2
            start[n]=$3
            end[n]=$4
            desc[n]=$5

            regions[1]="inside"
            win[1]=0
            regions[2]="500kb"
            win[2]=500000
            regions[3]="1Mb"
            win[3]=1000000
            regions[4]="2Mb"
            win[4]=2000000
        }

        best_file=OUT "/" NAME "_Campbell_best_pvalues.tsv"
        bonf_file=OUT "/" NAME "_Campbell_Bonferroni_hits.tsv"

        print "Gene","Gene_chr","Gene_start","Gene_end","Gene_description","Region","Window_bp","Best_SNP","Best_pos","Best_p_wald","Distance_bp","Direction","Bonferroni_threshold","Is_Bonferroni" > best_file

        print "Gene","Gene_chr","Gene_start","Gene_end","Gene_description","Region","Window_bp","SNP","SNP_pos","p_wald","Distance_bp","Direction","Bonferroni_threshold" > bonf_file
    }

    NR==1{next}

    {
        snp_chr=$1
        snp=$2
        pos=$3
        p=$13

        if(p=="nan" || p=="NA" || p=="" || p<=0 || p>1) next

        for(i=1;i<=n;i++){

            if(snp_chr != chr[i]) continue

            if(pos < start[i]){
                dist=start[i]-pos
                direction="upstream"
            } else if(pos > end[i]){
                dist=pos-end[i]
                direction="downstream"
            } else {
                dist=0
                direction="inside"
            }

            for(r=1;r<=4;r++){

                if(regions[r]=="inside" && dist != 0) continue
                if(regions[r]!="inside" && dist > win[r]) continue

                key=i SUBSEP r

                if(!(key in best) || p < best[key]){
                    best[key]=p
                    bestsnp[key]=snp
                    bestpos[key]=pos
                    bestdist[key]=dist
                    bestdir[key]=direction
                }

                if(p < BONF){
                    print gene[i],chr[i],start[i],end[i],desc[i],regions[r],win[r],snp,pos,p,dist,direction,BONF >> bonf_file
                }
            }
        }
    }

    END{
        for(i=1;i<=n;i++){
            for(r=1;r<=4;r++){
                key=i SUBSEP r

                if(key in best){
                    isbonf = (best[key] < BONF ? "YES" : "NO")
                    print gene[i],chr[i],start[i],end[i],desc[i],regions[r],win[r],bestsnp[key],bestpos[key],best[key],bestdist[key],bestdir[key],BONF,isbonf >> best_file
                } else {
                    print gene[i],chr[i],start[i],end[i],desc[i],regions[r],win[r],"NA","NA","NA","NA","NA",BONF,"NO" >> best_file
                }
            }
        }
    }
    ' "$GWAS"

    {
        echo "Analysis: $NAME"
        echo "GWAS file: $GWAS"
        echo "Output directory: $OUTDIR"
        echo "N tests: $NTESTS"
        echo "Bonferroni threshold: $BONF"
        echo
        echo "Best p-values:"
        cat "${OUTDIR}/${NAME}_Campbell_best_pvalues.tsv"
        echo
        echo "Number of Bonferroni hits near Campbell genes:"
        awk 'NR>1{n++} END{print n+0}' "${OUTDIR}/${NAME}_Campbell_Bonferroni_hits.tsv"
        echo
        echo "Bonferroni hits:"
        cat "${OUTDIR}/${NAME}_Campbell_Bonferroni_hits.tsv"
    } > "${OUTDIR}/${NAME}_Campbell_summary.txt"

    echo "Done: $NAME"
}

run_analysis \
"GEMMA_PC3_all" \
"${BASE}/result/GEMMA_PC3/output/tusk_GEMMA_LMM_PC3.assoc.txt" \
"${BASE}/tables_GEMMA_PC3/Campbell_candidates_v2"

run_analysis \
"GEMMA_male_PC3" \
"${BASE}/result/GEMMA_male/output/tusk_GEMMA_LMM_male_PC3.assoc.txt" \
"${BASE}/tables_GEMMA_male/Campbell_candidates_v2"

echo "All Campbell v2 analyses finished."
