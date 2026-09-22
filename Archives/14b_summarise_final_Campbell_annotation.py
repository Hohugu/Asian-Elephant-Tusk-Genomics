#!/usr/bin/env python3

from pathlib import Path
import csv
import re

BASE = Path("/scratch/project_2000886/Hoedric/GWAS_2025")
CAMP = BASE / "Genetics_Analysis" / "Campbell_exact_reproduction"
OUT = CAMP / "final_candidate_regions_annotation"

inp = OUT / "final_Campbell_candidate_regions_annotation.tsv"
out = OUT / "final_Campbell_candidate_regions_annotation_README_ready.tsv"

def extract_field(txt, key):
    m = re.search(rf"{key}=([^|;]+)", txt)
    return m.group(1) if m else "NA"

def simplify_overlap(txt):
    if txt == "no_overlapping_gene_like_feature":
        return "no direct overlap"

    genes = []
    for item in txt.split(";"):
        item = item.strip()
        if not item:
            continue
        ftype = item.split("|", 1)[0]
        name = extract_field(item, "Name")
        product = extract_field(item, "product")
        biotype = extract_field(item, "biotype")

        if ftype == "gene":
            genes.append((name, product, biotype))

    if genes:
        # unique while preserving order
        seen = set()
        out_items = []
        for name, product, biotype in genes:
            key = (name, product, biotype)
            if key not in seen:
                seen.add(key)
                if product != "NA":
                    out_items.append(f"{name} ({product}; {biotype})")
                else:
                    out_items.append(f"{name} ({biotype})")
        return "; ".join(out_items)

    return "overlaps transcript/exon/CDS but no gene feature parsed"

def simplify_nearest(txt):
    if txt.startswith("no_gene"):
        return "NA", "NA", "NA", "NA"

    ftype = txt.split("|", 1)[0]
    name = extract_field(txt, "Name")
    product = extract_field(txt, "product")
    biotype = extract_field(txt, "biotype")

    return ftype, name, product, biotype

with inp.open() as fh, out.open("w", newline="") as oh:
    reader = csv.DictReader(fh, delimiter="\t")

    fieldnames = [
        "region_id",
        "region",
        "signal",
        "threshold",
        "n_snps",
        "direct_overlap_summary",
        "nearest_gene_like_distance_bp",
        "nearest_gene_like_type",
        "nearest_gene_like_name",
        "nearest_gene_like_product",
        "nearest_gene_like_biotype",
        "README_annotation"
    ]

    writer = csv.DictWriter(oh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for row in reader:
        ftype, name, product, biotype = simplify_nearest(row["nearest_gene_like_feature"])
        overlap = simplify_overlap(row["overlapping_features"])

        region = f'{row["CHROM_CM"]}:{row["start"]}-{row["end"]}'

        if overlap != "no direct overlap":
            readme_annotation = overlap
        else:
            readme_annotation = (
                f"no direct overlap; nearest gene-like feature: "
                f"{name} ({ftype}, {biotype}), "
                f"{row['nearest_gene_like_distance_bp']} bp away"
            )

        writer.writerow({
            "region_id": row["region_id"],
            "region": region,
            "signal": row["signal"],
            "threshold": row["threshold"],
            "n_snps": row["n_snps"],
            "direct_overlap_summary": overlap,
            "nearest_gene_like_distance_bp": row["nearest_gene_like_distance_bp"],
            "nearest_gene_like_type": ftype,
            "nearest_gene_like_name": name,
            "nearest_gene_like_product": product,
            "nearest_gene_like_biotype": biotype,
            "README_annotation": readme_annotation
        })

print(out)
