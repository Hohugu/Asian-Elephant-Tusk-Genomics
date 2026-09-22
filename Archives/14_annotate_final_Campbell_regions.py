#!/usr/bin/env python3

from pathlib import Path
import csv
import re

BASE = Path("/scratch/project_2000886/Hoedric/GWAS_2025")
CAMP = BASE / "Genetics_Analysis" / "Campbell_exact_reproduction"

GFF = Path("/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff")
CM_NC = BASE / "GWAS" / "tables_GEMMA_male" / "CM_to_NC.map"

OUT = CAMP / "final_candidate_regions_annotation"
OUT.mkdir(parents=True, exist_ok=True)

regions = [
    {
        "region_id": "HetDev_CM044025_22Mb",
        "signal": "FST+DXY+HetDev_TX",
        "threshold": "5pct_projected/strictHet",
        "CHROM_CM": "CM044025.1",
        "start": 22575792,
        "end": 22576805,
        "n_snps": "20_projected_2_strictHet",
        "comment": "HetDev-overlap sensitivity signal",
    },
    {
        "region_id": "FDM_CM044025_143Mb",
        "signal": "FST+DXY+muLD",
        "threshold": "1pct",
        "CHROM_CM": "CM044025.1",
        "start": 143590001,
        "end": 143600000,
        "n_snps": "35",
        "comment": "Main strict Campbell-like signal",
    },
    {
        "region_id": "FDM_CM044047_34Mb",
        "signal": "FST+DXY+muLD",
        "threshold": "1pct",
        "CHROM_CM": "CM044047.1",
        "start": 34830001,
        "end": 34850000,
        "n_snps": "30",
        "comment": "Merged adjacent bins",
    },
    {
        "region_id": "FDM_CM044022_26Mb",
        "signal": "FST+DXY+muLD",
        "threshold": "1pct",
        "CHROM_CM": "CM044022.1",
        "start": 26030001,
        "end": 26040000,
        "n_snps": "18",
        "comment": "Main strict Campbell-like signal",
    },
    {
        "region_id": "FDM_CM044020_70Mb",
        "signal": "FST+DXY+muLD",
        "threshold": "1pct",
        "CHROM_CM": "CM044020.1",
        "start": 70690001,
        "end": 70700000,
        "n_snps": "5",
        "comment": "Main strict Campbell-like signal",
    },
    {
        "region_id": "FDM_CM044023_37Mb",
        "signal": "FST+DXY+muLD",
        "threshold": "1pct",
        "CHROM_CM": "CM044023.1",
        "start": 37750001,
        "end": 37760000,
        "n_snps": "4",
        "comment": "Main strict Campbell-like signal",
    },
]

def parse_attrs(attr):
    out = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out

def clean_value(x):
    if x is None or x == "":
        return "NA"
    return x.replace("\t", " ").replace("\n", " ")

# Load CM -> NC mapping
cm_to_nc = {}
with CM_NC.open() as f:
    for line in f:
        if not line.strip():
            continue
        cm, nc = line.rstrip("\n").split("\t")[:2]
        cm_to_nc[cm] = nc

for r in regions:
    r["CHROM_NC"] = cm_to_nc.get(r["CHROM_CM"], "NA")
    r["overlapping_features"] = []
    r["nearest_gene_like_distance_bp"] = None
    r["nearest_gene_like_feature"] = None

# Write region input table
region_file = OUT / "final_Campbell_candidate_regions_with_mapping.tsv"
with region_file.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "region_id", "signal", "threshold", "CHROM_CM", "CHROM_NC",
            "start", "end", "n_snps", "comment"
        ],
        delimiter="\t",
    )
    writer.writeheader()
    for r in regions:
        writer.writerow({k: r[k] for k in writer.fieldnames})

# Check scaffold presence in GFF
seqid_counts = {r["CHROM_NC"]: 0 for r in regions}

# We use these as overlapping annotation features.
overlap_types = {
    "gene", "pseudogene", "lnc_RNA", "mRNA", "transcript",
    "exon", "CDS", "ncRNA", "rRNA", "tRNA"
}

# We use these for nearest gene-like feature.
nearest_types = {"gene", "pseudogene", "lnc_RNA", "ncRNA"}

with GFF.open() as f:
    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        seqid, source, ftype, start, end, score, strand, phase, attr = parts

        try:
            start = int(start)
            end = int(end)
        except ValueError:
            continue

        if seqid in seqid_counts:
            seqid_counts[seqid] += 1

        if ftype == "region":
            continue

        attrs = parse_attrs(attr)

        gene_id = clean_value(attrs.get("ID", attrs.get("Dbxref", "NA")))
        gene_name = clean_value(attrs.get("Name", attrs.get("gene", "NA")))
        product = clean_value(attrs.get("product", attrs.get("description", "NA")))
        biotype = clean_value(attrs.get("gene_biotype", attrs.get("gbkey", "NA")))

        feature_summary = (
            f"{ftype}|ID={gene_id}|Name={gene_name}|"
            f"product={product}|biotype={biotype}|"
            f"{seqid}:{start}-{end}|strand={strand}"
        )

        for r in regions:
            if seqid != r["CHROM_NC"]:
                continue

            rstart = int(r["start"])
            rend = int(r["end"])

            # Overlap
            if ftype in overlap_types and start <= rend and end >= rstart:
                r["overlapping_features"].append(feature_summary)

            # Nearest gene-like
            if ftype in nearest_types:
                if end < rstart:
                    dist = rstart - end
                elif start > rend:
                    dist = start - rend
                else:
                    dist = 0

                if (
                    r["nearest_gene_like_distance_bp"] is None
                    or dist < r["nearest_gene_like_distance_bp"]
                ):
                    r["nearest_gene_like_distance_bp"] = dist
                    r["nearest_gene_like_feature"] = feature_summary

# Write scaffold presence
presence_file = OUT / "scaffold_presence_in_GFF_after_CM_to_NC_mapping.tsv"
with presence_file.open("w") as f:
    f.write("CHROM_NC\tn_features_in_GFF\n")
    for seqid, count in sorted(seqid_counts.items()):
        f.write(f"{seqid}\t{count}\n")

# Write final annotation table
annotation_file = OUT / "final_Campbell_candidate_regions_annotation.tsv"

with annotation_file.open("w", newline="") as f:
    fieldnames = [
        "region_id",
        "signal",
        "threshold",
        "CHROM_CM",
        "CHROM_NC",
        "start",
        "end",
        "n_snps",
        "overlapping_features",
        "nearest_gene_like_distance_bp",
        "nearest_gene_like_feature",
        "comment",
    ]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for r in regions:
        overlapping = r["overlapping_features"]
        if overlapping:
            overlapping_text = "; ".join(overlapping)
        else:
            overlapping_text = "no_overlapping_gene_like_feature"

        nearest_dist = r["nearest_gene_like_distance_bp"]
        nearest_feature = r["nearest_gene_like_feature"]

        if nearest_dist is None:
            nearest_dist = "NA"
            nearest_feature = "no_gene_like_feature_on_scaffold"

        writer.writerow({
            "region_id": r["region_id"],
            "signal": r["signal"],
            "threshold": r["threshold"],
            "CHROM_CM": r["CHROM_CM"],
            "CHROM_NC": r["CHROM_NC"],
            "start": r["start"],
            "end": r["end"],
            "n_snps": r["n_snps"],
            "overlapping_features": overlapping_text,
            "nearest_gene_like_distance_bp": nearest_dist,
            "nearest_gene_like_feature": nearest_feature,
            "comment": r["comment"],
        })

print("Saved:")
print(region_file)
print(presence_file)
print(annotation_file)
