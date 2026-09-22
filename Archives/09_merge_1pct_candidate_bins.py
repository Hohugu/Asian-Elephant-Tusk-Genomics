#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import Counter

OUT = Path("/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction/snp_level_campbell")

infile = OUT / "Campbell_ge3_1pct_candidate_bins_summary.tsv"
outfile = OUT / "Campbell_ge3_1pct_candidate_regions_merged.tsv"

rows = []

with open(infile, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        row["bin_start"] = int(row["bin_start"])
        row["bin_end"] = int(row["bin_end"])
        row["n_snps"] = int(row["n_snps"])
        rows.append(row)

rows.sort(key=lambda r: (r["CHROM"], r["bin_start"], r["bin_end"]))

merged = []

for row in rows:
    if (
        not merged
        or row["CHROM"] != merged[-1]["CHROM"]
        or row["bin_start"] > merged[-1]["region_end"] + 1
    ):
        merged.append({
            "CHROM": row["CHROM"],
            "region_start": row["bin_start"],
            "region_end": row["bin_end"],
            "n_bins": 1,
            "n_snps": row["n_snps"],
            "bins": [row["bin_id"]],
            "patterns": Counter(),
            "classes": Counter(),
            "overlap_genes": set(),
            "nearest_genes": set(),
            "nearest_gene_descriptions": set(),
        })
    else:
        merged[-1]["region_end"] = max(merged[-1]["region_end"], row["bin_end"])
        merged[-1]["n_bins"] += 1
        merged[-1]["n_snps"] += row["n_snps"]
        merged[-1]["bins"].append(row["bin_id"])

    m = merged[-1]

    for item in row["patterns"].split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            m["patterns"][k] += int(v)

    for item in row["classes"].split(";"):
        if ":" in item:
            k, v = item.split(":", 1)
            m["classes"][k] += int(v)

    for x in row["overlap_genes"].split(","):
        if x:
            m["overlap_genes"].add(x)

    for x in row["nearest_genes"].split(","):
        if x:
            m["nearest_genes"].add(x)

    for x in row["nearest_gene_descriptions"].split(","):
        if x:
            m["nearest_gene_descriptions"].add(x)

def fmt_counter(c):
    return ";".join(f"{k}:{v}" for k, v in c.most_common())

def fmt_set(s):
    return ",".join(sorted(s))

merged.sort(key=lambda r: r["n_snps"], reverse=True)

with open(outfile, "w", newline="") as out:
    fieldnames = [
        "region_id",
        "CHROM",
        "region_start",
        "region_end",
        "region_size_bp",
        "n_bins",
        "n_snps",
        "patterns",
        "classes",
        "overlap_genes",
        "nearest_genes",
        "nearest_gene_descriptions",
        "bins"
    ]

    writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for i, r in enumerate(merged, start=1):
        writer.writerow({
            "region_id": f"Campbell_1pct_region_{i}",
            "CHROM": r["CHROM"],
            "region_start": r["region_start"],
            "region_end": r["region_end"],
            "region_size_bp": r["region_end"] - r["region_start"] + 1,
            "n_bins": r["n_bins"],
            "n_snps": r["n_snps"],
            "patterns": fmt_counter(r["patterns"]),
            "classes": fmt_counter(r["classes"]),
            "overlap_genes": fmt_set(r["overlap_genes"]),
            "nearest_genes": fmt_set(r["nearest_genes"]),
            "nearest_gene_descriptions": fmt_set(r["nearest_gene_descriptions"]),
            "bins": ",".join(r["bins"])
        })

print(f"Wrote {outfile}")
