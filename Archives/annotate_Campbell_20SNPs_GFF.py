#!/usr/bin/env python3

from pathlib import Path
import gzip
import re

BASE = Path("/scratch/project_2000886/Hoedric/GWAS_2025")
CAMP = BASE / "Genetics_Analysis/Campbell_exact_reproduction"
AUDIT = CAMP / "HetDev_audit"

SNP_FILE = AUDIT / "Campbell_5pct_exact_FST_DXY_HetDev_TX_20SNPs.tsv"
GFF = Path("/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff")

OUT = AUDIT / "Campbell_5pct_exact_FST_DXY_HetDev_TX_20SNPs_GFF_annotation.tsv"

def parse_attrs(s):
    d = {}
    for item in s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def feature_name(attrs):
    for k in ["Name", "gene", "gene_name", "product", "ID", "Parent"]:
        if k in attrs:
            return attrs[k]
    return ""

# Read SNPs
snps = []
with open(SNP_FILE) as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {x:i for i,x in enumerate(header)}
    for line in f:
        row = line.rstrip("\n").split("\t")
        chrom = row[idx["CHROM"]]
        pos = int(float(row[idx["POS"]]))
        snps.append({
            "snp_id": row[idx["snp_id"]],
            "CHROM": chrom,
            "POS": pos,
            "REF": row[idx["REF"]],
            "ALT": row[idx["ALT"]],
            "pattern": row[idx["pattern"]],
            "DXY_value": row[idx["DXY_value"]],
            "HetDev_TX_value": row[idx["HetDev_TX_value"]],
            "TX_p_alt": row[idx["TX_p_alt"]],
            "TT_alt_sum": row[idx["TT_alt_sum"]],
            "Campbell_Het_candidate": row[idx["Campbell_Het_candidate"]],
        })

snps_by_chr = {}
for s in snps:
    snps_by_chr.setdefault(s["CHROM"], []).append(s)

# Read GFF features
genes_by_chr = {}
features_by_chr = {}

opener = gzip.open if str(GFF).endswith(".gz") else open
with opener(GFF, "rt") as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, source, ftype, start, end, score, strand, phase, attrs = parts
        if chrom not in snps_by_chr:
            continue
        start = int(start)
        end = int(end)
        attrd = parse_attrs(attrs)
        name = feature_name(attrd)
        rec = {
            "CHROM": chrom,
            "type": ftype,
            "start": start,
            "end": end,
            "strand": strand,
            "name": name,
            "attrs": attrs
        }
        features_by_chr.setdefault(chrom, []).append(rec)
        if ftype in ["gene", "pseudogene", "lnc_RNA", "transcript", "mRNA"]:
            genes_by_chr.setdefault(chrom, []).append(rec)

# Annotate
out_rows = []

for s in snps:
    chrom = s["CHROM"]
    pos = s["POS"]

    overlapping = []
    for feat in features_by_chr.get(chrom, []):
        if feat["start"] <= pos <= feat["end"]:
            overlapping.append(feat)

    overlapping_genes = [
        f for f in overlapping
        if f["type"] in ["gene", "pseudogene", "lnc_RNA", "transcript", "mRNA"]
    ]

    # classify
    if overlapping:
        feature_types = ",".join(sorted(set(f["type"] for f in overlapping)))
        feature_names = ",".join(sorted(set(f["name"] for f in overlapping if f["name"])))
    else:
        feature_types = "intergenic"
        feature_names = ""

    # nearest gene-like feature
    nearest = None
    nearest_dist = None
    for g in genes_by_chr.get(chrom, []):
        if g["start"] <= pos <= g["end"]:
            dist = 0
        elif pos < g["start"]:
            dist = g["start"] - pos
        else:
            dist = pos - g["end"]

        if nearest_dist is None or dist < nearest_dist:
            nearest_dist = dist
            nearest = g

    if nearest is None:
        nearest_name = ""
        nearest_type = ""
        nearest_start = ""
        nearest_end = ""
        nearest_strand = ""
        nearest_dist = ""
        relation = "no_gene_on_scaffold"
    else:
        nearest_name = nearest["name"]
        nearest_type = nearest["type"]
        nearest_start = nearest["start"]
        nearest_end = nearest["end"]
        nearest_strand = nearest["strand"]

        if nearest_dist == 0:
            relation = "inside_gene_like_feature"
        elif nearest_dist <= 2000:
            relation = "within_2kb_of_gene_like_feature"
        elif nearest_dist <= 100000:
            relation = "within_100kb_of_gene_like_feature"
        else:
            relation = "intergenic_distal"

    out_rows.append({
        **s,
        "overlap_feature_types": feature_types,
        "overlap_feature_names": feature_names,
        "nearest_gene_like_name": nearest_name,
        "nearest_gene_like_type": nearest_type,
        "nearest_gene_like_start": nearest_start,
        "nearest_gene_like_end": nearest_end,
        "nearest_gene_like_strand": nearest_strand,
        "distance_to_nearest_gene_like": nearest_dist,
        "relation_to_nearest_gene_like": relation
    })

cols = [
    "snp_id","CHROM","POS","REF","ALT","pattern",
    "DXY_value","HetDev_TX_value","TX_p_alt","TT_alt_sum","Campbell_Het_candidate",
    "overlap_feature_types","overlap_feature_names",
    "nearest_gene_like_name","nearest_gene_like_type",
    "nearest_gene_like_start","nearest_gene_like_end","nearest_gene_like_strand",
    "distance_to_nearest_gene_like","relation_to_nearest_gene_like"
]

with open(OUT, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in out_rows:
        out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

print("Wrote:", OUT)
print("n_snps:", len(out_rows))
