#!/bin/bash
set -euo pipefail

BASE=/scratch/project_2000886/Hoedric/GWAS_2025
CAMP=$BASE/Genetics_Analysis/Campbell_exact_reproduction
GFF=/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff

OUT=$CAMP/final_candidate_regions_annotation
mkdir -p "$OUT"

REGIONS=$OUT/final_Campbell_candidate_regions.tsv

cat > "$REGIONS" <<'EOR'
region_id	signal	threshold	CHROM	start	end	n_snps	comment
HetDev_CM044025_22Mb	FST+DXY+HetDev_TX	5pct_projected/strictHet	CM044025.1	22575792	22576805	20_projected_2_strictHet	HetDev-overlap sensitivity signal
FDM_CM044025_143Mb	FST+DXY+muLD	1pct	CM044025.1	143590001	143600000	35	Main strict Campbell-like signal
FDM_CM044047_34Mb	FST+DXY+muLD	1pct	CM044047.1	34830001	34850000	30	Merged adjacent bins
FDM_CM044022_26Mb	FST+DXY+muLD	1pct	CM044022.1	26030001	26040000	18	Main strict Campbell-like signal
FDM_CM044020_70Mb	FST+DXY+muLD	1pct	CM044020.1	70690001	70700000	5	Main strict Campbell-like signal
FDM_CM044023_37Mb	FST+DXY+muLD	1pct	CM044023.1	37750001	37760000	4	Main strict Campbell-like signal
EOR

echo "=== Check if scaffolds exist in GFF ==="
cut -f4 "$REGIONS" | tail -n +2 | sort -u | while read chr; do
  n=$(grep -v '^#' "$GFF" | awk -v c="$chr" '$1==c{n++} END{print n+0}')
  echo -e "$chr\t$n"
done > "$OUT/scaffold_presence_in_GFF.tsv"

cat "$OUT/scaffold_presence_in_GFF.tsv"

echo
echo "=== Example sequence names in GFF ==="
grep -v '^#' "$GFF" | awk -F'\t' '{print $1}' | sort -u | head -20 > "$OUT/example_GFF_seqids.txt"
cat "$OUT/example_GFF_seqids.txt"

echo
echo "=== Annotating overlapping gene-like features and nearest genes ==="

awk -F'\t' -v OFS='\t' '
BEGIN {
  first = 1
  while ((getline line < ARGV[1]) > 0) {
    if (first == 1) {
      first = 0
      continue
    }

    split(line, a, "\t")
    rid=a[1]; sig=a[2]; thr=a[3]; chr=a[4]; start=a[5]; end=a[6]; nsnp=a[7]; comment=a[8]

    region_chr[rid]=chr
    region_start[rid]=start
    region_end[rid]=end
    region_sig[rid]=sig
    region_thr[rid]=thr
    region_nsnp[rid]=nsnp
    region_comment[rid]=comment
    region_ids[++n_regions]=rid
  }

  ARGV[1] = ""
}
$0 !~ /^#/ {
  chr=$1
  type=$3
  start=$4
  end=$5
  attr=$9

  if (type ~ /^(gene|mRNA|transcript|lnc_RNA|pseudogene|exon|CDS)$/) {

    gene_id="NA"
    gene_name="NA"
    product="NA"

    if (match(attr, /ID=([^;]+)/)) {
      gene_id=substr(attr, RSTART+3, RLENGTH-3)
    }

    if (match(attr, /Name=([^;]+)/)) {
      gene_name=substr(attr, RSTART+5, RLENGTH-5)
    }

    if (match(attr, /gene=([^;]+)/)) {
      gene_name=substr(attr, RSTART+5, RLENGTH-5)
    }

    if (match(attr, /product=([^;]+)/)) {
      product=substr(attr, RSTART+8, RLENGTH-8)
    }

    for (i=1; i<=n_regions; i++) {
      rid=region_ids[i]
      rchr=region_chr[rid]
      rstart=region_start[rid]
      rend=region_end[rid]

      if (chr == rchr) {

        if (start <= rend && end >= rstart) {
          overlap[rid] = overlap[rid] type ":" gene_id ":" gene_name ":" product ":" start "-" end ";"
        }

        if (type ~ /^(gene|lnc_RNA|pseudogene)$/) {
          if (end < rstart) {
            dist = rstart - end
          } else if (start > rend) {
            dist = start - rend
          } else {
            dist = 0
          }

          if (!(rid in bestdist) || dist < bestdist[rid]) {
            bestdist[rid] = dist
            bestfeature[rid] = type ":" gene_id ":" gene_name ":" product ":" start "-" end
          }
        }
      }
    }
  }
}
END {
  print "region_id","signal","threshold","CHROM","start","end","n_snps","overlapping_features","nearest_gene_like_distance_bp","nearest_gene_like_feature","comment"

  for (i=1; i<=n_regions; i++) {
    rid=region_ids[i]

    ov = overlap[rid]
    if (ov == "") ov = "no_overlapping_gene_like_feature_or_seqid_mismatch"

    if (rid in bestdist) {
      bd = bestdist[rid]
      bf = bestfeature[rid]
    } else {
      bd = "NA"
      bf = "no_gene_like_feature_on_scaffold_or_seqid_mismatch"
    }

    print rid, region_sig[rid], region_thr[rid], region_chr[rid], region_start[rid], region_end[rid], region_nsnp[rid], ov, bd, bf, region_comment[rid]
  }
}
' "$REGIONS" "$GFF" > "$OUT/final_Campbell_candidate_regions_annotation.tsv"

echo
echo "=== Output ==="
column -t -s $'\t' "$OUT/final_Campbell_candidate_regions_annotation.tsv"

echo
echo "Saved:"
echo "$OUT/final_Campbell_candidate_regions_annotation.tsv"
