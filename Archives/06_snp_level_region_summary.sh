#!/bin/bash
set -euo pipefail

export BASE=/scratch/project_2000886/Hoedric/GWAS_2025
export CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction
export OUT=$CAMP/snp_level_campbell

SUMMARY=$OUT/Campbell_SNP_level_priority_region_summary.tsv
PATTERNS=$OUT/Campbell_SNP_level_priority_region_patterns.tsv

printf "threshold\tregion\tCHROM\tstart\tend\tn_any\tn_FST\tn_DXY\tn_muLD\tn_HetDev_TX\tn_FST_DXY\tn_FST_muLD\tn_FST_HetDev_TX\tn_DXY_muLD\tn_DXY_HetDev_TX\tn_muLD_HetDev_TX\tn_FST_DXY_muLD\tn_FST_DXY_HetDev_TX\tn_FST_muLD_HetDev_TX\tn_DXY_muLD_HetDev_TX\tn_all4\n" > "$SUMMARY"

printf "threshold\tregion\tpattern\tcount\n" > "$PATTERNS"

summarize_region () {
  local threshold=$1
  local infile=$2
  local region=$3
  local chr=$4
  local start=$5
  local end=$6

  zcat "$infile" | awk -F'\t' \
    -v OFS='\t' \
    -v threshold="$threshold" \
    -v region="$region" \
    -v chr="$chr" \
    -v start="$start" \
    -v end="$end" \
    '
    NR>1 && $2==chr && $3>=start && $3<=end {
      F=($6=="True")
      D=($7=="True")
      M=($8=="True")
      H=($9=="True")

      n_any++
      n_FST += F
      n_DXY += D
      n_muLD += M
      n_Het += H

      n_FST_DXY += (F && D)
      n_FST_muLD += (F && M)
      n_FST_Het += (F && H)
      n_DXY_muLD += (D && M)
      n_DXY_Het += (D && H)
      n_muLD_Het += (M && H)

      n_FST_DXY_muLD += (F && D && M)
      n_FST_DXY_Het += (F && D && H)
      n_FST_muLD_Het += (F && M && H)
      n_DXY_muLD_Het += (D && M && H)
      n_all4 += (F && D && M && H)

      pattern[$11]++
    }
    END {
      print threshold, region, chr, start, end,
            n_any+0,
            n_FST+0, n_DXY+0, n_muLD+0, n_Het+0,
            n_FST_DXY+0, n_FST_muLD+0, n_FST_Het+0,
            n_DXY_muLD+0, n_DXY_Het+0, n_muLD_Het+0,
            n_FST_DXY_muLD+0,
            n_FST_DXY_Het+0,
            n_FST_muLD_Het+0,
            n_DXY_muLD_Het+0,
            n_all4+0
    }' >> "$SUMMARY"

  zcat "$infile" | awk -F'\t' \
    -v OFS='\t' \
    -v threshold="$threshold" \
    -v region="$region" \
    -v chr="$chr" \
    -v start="$start" \
    -v end="$end" \
    '
    NR>1 && $2==chr && $3>=start && $3<=end {
      pattern[$11]++
    }
    END {
      for (p in pattern) print threshold, region, p, pattern[p]
    }' >> "$PATTERNS"
}

FILE5=$OUT/Campbell_SNP_membership_5pct_Campbell.tsv.gz
FILE1=$OUT/Campbell_SNP_membership_1pct_Campbell_threshold.tsv.gz

for threshold in 5pct_Campbell 1pct_Campbell_threshold
do
  if [ "$threshold" = "5pct_Campbell" ]; then
    infile=$FILE5
  else
    infile=$FILE1
  fi

  summarize_region "$threshold" "$infile" "LOC126069858_gene_body" "CM044048.1" 14929369 15115399
  summarize_region "$threshold" "$infile" "LOC126069858_boundary_10kb_bin" "CM044048.1" 14920001 14930000
  summarize_region "$threshold" "$infile" "LOC126069858_pm100kb" "CM044048.1" 14829369 15215399

  summarize_region "$threshold" "$infile" "AMELX_gene_body" "CM044047.1" 168725166 168729418
  summarize_region "$threshold" "$infile" "AMELX_regional_window" "CM044047.1" 167500000 171000000
done

echo "Wrote:"
echo "$SUMMARY"
echo "$PATTERNS"
