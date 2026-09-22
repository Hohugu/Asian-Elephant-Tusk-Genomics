#!/usr/bin/env python3

import csv
from collections import defaultdict, Counter
from pathlib import Path

OUT = Path("/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/snp_level_campbell")

infile = OUT / "Campbell_ge3_SNPs_1pct_Campbell_threshold.annotated.tsv"
outfile = OUT / "Campbell_ge3_1pct_candidate_bins_summary.tsv"

bins = {}

with open(infile, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        bin_id = row["bin_id"]
        pos = int(row["POS"])

        if bin_id not in bins:
            bins[bin_id] = {
                "bin_id": bin_id,
                "CHROM": row["CHROM"],
                "bin_start": row["bin_start"],
                "bin_end": row["bin_end"],
                "n_snps": 0,
                "min_POS": pos,
                "max_POS": pos,
                "patterns": Counter(),
                "classes": Counter(),
                "overlap_genes": set(),
                "nearest_genes": set(),
                "nearest_gene_descriptions": set()
            }

        b = bins[bin_id]
        b["n_snps"] += 1
        b["min_POS"] = min(b["min_POS"], pos)
        b["max_POS"] = max(b["max_POS"], pos)

        b["patterns"][row["pattern"]] += 1
        b["classes"][row["genomic_class"]] += 1

        if row["overlap_gene"]:
            b["overlap_genes"].add(row["overlap_gene"])

        if row["nearest_gene"]:
            b["nearest_genes"].add(row["nearest_gene"])

        if row["nearest_gene_description"]:
            b["nearest_gene_descriptions"].add(row["nearest_gene_description"])

def format_counter(c):
    return ";".join(f"{k}:{v}" for k, v in c.most_common())

def format_set(s):
    return ",".join(sorted(s))

rows = sorted(
    bins.values(),
    key=lambda x: x["n_snps"],
    reverse=True
)

with open(outfile, "w", newline="") as out:
    fieldnames = [
        "bin_id",
        "CHROM",
        "bin_start",
        "bin_end",
        "n_snps",
        "min_POS",
        "max_POS",
        "patterns",
        "classes",
        "overlap_genes",
        "nearest_genes",
        "nearest_gene_descriptions"
    ]

    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for b in rows:
        writer.writerow({
            "bin_id": b["bin_id"],
            "CHROM": b["CHROM"],
            "bin_start": b["bin_start"],
            "bin_end": b["bin_end"],
            "n_snps": b["n_snps"],
            "min_POS": b["min_POS"],
            "max_POS": b["max_POS"],
            "patterns": format_counter(b["patterns"]),
            "classes": format_counter(b["classes"]),
            "overlap_genes": format_set(b["overlap_genes"]),
            "nearest_genes": format_set(b["nearest_genes"]),
            "nearest_gene_descriptions": format_set(b["nearest_gene_descriptions"])
        })

print(f"Wrote {outfile}")
