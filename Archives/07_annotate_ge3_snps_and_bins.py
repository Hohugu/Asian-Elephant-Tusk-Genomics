#!/usr/bin/env python3

import csv
import gzip
from pathlib import Path
from collections import defaultdict, Counter

CAMP = Path("/scratch/project_2000886/Hoedric/GWAS_2025/Genetics_Analysis/Campbell_exact_reproduction")
OUT = CAMP / "snp_level_campbell"
GFF = Path("/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff")

PROMOTER_BP = 2000

INPUTS = [
    ("5pct_Campbell", OUT / "Campbell_SNPs_ge3analyses_5pct.tsv"),
    ("1pct_Campbell_threshold", OUT / "Campbell_SNPs_ge3analyses_1pct.tsv"),
]

def nc_to_cm(seqid):
    """
    Convert NC_064847.1 -> CM044048.1 for primary chromosomes/scaffolds.
    Keeps non-matching contigs unchanged.
    """
    if seqid.startswith("NC_064") and seqid.endswith(".1"):
        try:
            n = int(seqid.replace("NC_", "").replace(".1", ""))
            cm = n - 20799
            if 44000 <= cm <= 45000:
                return f"CM{cm:06d}.1"
        except ValueError:
            pass
    return seqid

def parse_attrs(attr):
    d = {}
    for item in attr.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def load_gff_features(gff):
    features = {
        "gene": defaultdict(list),
        "exon": defaultdict(list),
        "CDS": defaultdict(list),
        "promoter_2kb": defaultdict(list),
    }

    with open_text(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attr = parts
            chrom = nc_to_cm(seqid)
            start = int(start)
            end = int(end)
            attrs = parse_attrs(attr)

            gene_name = attrs.get("gene", attrs.get("Name", attrs.get("ID", "")))
            gene_id = attrs.get("ID", "")
            description = attrs.get("description", "")
            biotype = attrs.get("gene_biotype", "")

            record = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "gene": gene_name,
                "gene_id": gene_id,
                "description": description,
                "biotype": biotype,
                "feature": ftype,
            }

            if ftype == "gene":
                features["gene"][chrom].append(record)

                tss = start if strand != "-" else end
                p_start = max(1, tss - PROMOTER_BP)
                p_end = tss + PROMOTER_BP

                prom = dict(record)
                prom["start"] = p_start
                prom["end"] = p_end
                prom["feature"] = "promoter_2kb"
                features["promoter_2kb"][chrom].append(prom)

            elif ftype == "exon":
                features["exon"][chrom].append(record)

            elif ftype == "CDS":
                features["CDS"][chrom].append(record)

    for ftype in features:
        for chrom in features[ftype]:
            features[ftype][chrom].sort(key=lambda x: (x["start"], x["end"]))

    return features

def overlap_records(records, pos):
    hits = []
    for r in records:
        if r["start"] > pos:
            break
        if r["end"] >= pos:
            hits.append(r)
    return hits

def nearest_gene(gene_records, pos):
    if not gene_records:
        return "", "", ""

    best = None
    best_dist = None

    for g in gene_records:
        if g["start"] <= pos <= g["end"]:
            return g["gene"], "0", g["description"]

        if pos < g["start"]:
            dist = g["start"] - pos
        else:
            dist = pos - g["end"]

        if best_dist is None or dist < best_dist:
            best = g
            best_dist = dist

    return best["gene"], str(best_dist), best["description"]

def unique_join(values):
    vals = []
    seen = set()
    for v in values:
        if v and v not in seen:
            vals.append(v)
            seen.add(v)
    return ",".join(vals)

def annotate_snp(row, features):
    chrom = row["CHROM"]
    pos = int(row["POS"])

    gene_hits = overlap_records(features["gene"].get(chrom, []), pos)
    exon_hits = overlap_records(features["exon"].get(chrom, []), pos)
    cds_hits = overlap_records(features["CDS"].get(chrom, []), pos)
    prom_hits = overlap_records(features["promoter_2kb"].get(chrom, []), pos)

    in_gene = len(gene_hits) > 0
    in_exon = len(exon_hits) > 0
    in_cds = len(cds_hits) > 0
    in_prom = len(prom_hits) > 0

    if in_cds:
        genomic_class = "CDS"
    elif in_exon:
        genomic_class = "exon_non_CDS"
    elif in_prom:
        genomic_class = "promoter_2kb"
    elif in_gene:
        genomic_class = "intragenic_non_exonic"
    else:
        genomic_class = "intergenic"

    overlap_genes = gene_hits if gene_hits else prom_hits
    overlap_gene_names = unique_join([g["gene"] for g in overlap_genes])
    overlap_gene_desc = unique_join([g["description"] for g in overlap_genes])

    nearest_name, nearest_dist, nearest_desc = nearest_gene(features["gene"].get(chrom, []), pos)

    bin_start = ((pos - 1) // 10000) * 10000 + 1
    bin_end = bin_start + 9999
    bin_id = f"{chrom}:{bin_start}:{bin_end}"

    out = dict(row)
    out.update({
        "bin_id": bin_id,
        "bin_start": bin_start,
        "bin_end": bin_end,
        "genomic_class": genomic_class,
        "in_gene": in_gene,
        "in_exon": in_exon,
        "in_CDS": in_cds,
        "in_promoter_2kb": in_prom,
        "overlap_gene": overlap_gene_names,
        "overlap_gene_description": overlap_gene_desc,
        "nearest_gene": nearest_name,
        "nearest_gene_distance_bp": nearest_dist,
        "nearest_gene_description": nearest_desc,
    })

    return out

print("Loading GFF features", flush=True)
features = load_gff_features(GFF)

all_summary = []
all_bin_rows = []

for threshold, infile in INPUTS:
    print(f"Annotating {threshold}: {infile}", flush=True)

    annotated_rows = []

    with open(infile, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            annotated_rows.append(annotate_snp(row, features))

    snp_out = OUT / f"Campbell_ge3_SNPs_{threshold}.annotated.tsv"

    fieldnames = list(annotated_rows[0].keys()) if annotated_rows else []

    with open(snp_out, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(annotated_rows)

    class_counts = Counter(r["genomic_class"] for r in annotated_rows)
    pattern_counts = Counter(r["pattern"] for r in annotated_rows)

    for k, v in class_counts.items():
        all_summary.append({
            "threshold": threshold,
            "summary_type": "genomic_class",
            "category": k,
            "count": v
        })

    for k, v in pattern_counts.items():
        all_summary.append({
            "threshold": threshold,
            "summary_type": "pattern",
            "category": k,
            "count": v
        })

    bins = defaultdict(list)
    for r in annotated_rows:
        bins[r["bin_id"]].append(r)

    for bin_id, rows in bins.items():
        first = rows[0]
        patterns = Counter(r["pattern"] for r in rows)
        classes = Counter(r["genomic_class"] for r in rows)

        all_bin_rows.append({
            "threshold": threshold,
            "bin_id": bin_id,
            "CHROM": first["CHROM"],
            "bin_start": first["bin_start"],
            "bin_end": first["bin_end"],
            "n_ge3_snps": len(rows),
            "patterns": ";".join(f"{k}:{v}" for k, v in patterns.most_common()),
            "classes": ";".join(f"{k}:{v}" for k, v in classes.most_common()),
            "overlap_genes": unique_join([r["overlap_gene"] for r in rows]),
            "nearest_genes": unique_join([r["nearest_gene"] for r in rows]),
        })

    print(f"Wrote {snp_out}", flush=True)

summary_out = OUT / "Campbell_ge3_annotation_summary.tsv"
with open(summary_out, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t",
                            fieldnames=["threshold", "summary_type", "category", "count"])
    writer.writeheader()
    writer.writerows(all_summary)

bins_out = OUT / "Campbell_ge3_bins_10kb.annotated.tsv"
with open(bins_out, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t",
                            fieldnames=[
                                "threshold", "bin_id", "CHROM", "bin_start", "bin_end",
                                "n_ge3_snps", "patterns", "classes",
                                "overlap_genes", "nearest_genes"
                            ])
    writer.writeheader()
    writer.writerows(all_bin_rows)

print(f"Wrote {summary_out}", flush=True)
print(f"Wrote {bins_out}", flush=True)
