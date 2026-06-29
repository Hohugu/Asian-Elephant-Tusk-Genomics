#!/usr/bin/env python3

import os
import re
import subprocess
import tempfile
from collections import defaultdict

BASE = "/scratch/project_2000886/Hoedric/GWAS_2025"

GFF = "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

FASTA_CANDIDATES = [
    "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna",
    "/scratch/project_2000886/Elisa/reseq/GCA_024166365.1_mEleMax1_primary_haplotype_genomic.fa",
    "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fa",
]

MAP = f"{BASE}/GWAS/tables_GEMMA_male/CM_to_NC.map"

RAW_VCF = f"{BASE}/results/genotyped_vcf/all_contigs.raw.vcf.gz"

RARE_ANNOT = (
    f"{BASE}/Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/"
    "TT_TX_rawVCF_rare_variants_regulatory_priority_annotation.tsv"
)

OUT = (
    f"{BASE}/Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/"
    "TT_TX_rawVCF_rare_CDS_variant_effects.tsv"
)

SUMMARY = (
    f"{BASE}/Genetics_Analysis/Functional_annotation/rare_variants_rawVCF/"
    "TT_TX_rawVCF_rare_CDS_variant_effects_summary.tsv"
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def run(cmd):
    return subprocess.check_output(cmd, text=True)


def revcomp(seq):
    return seq.translate(COMP)[::-1].upper()


def parse_attr(attr, key):
    match = re.search(rf"{key}=([^;]+)", attr)
    if match:
        return match.group(1)
    return ""


def decode_gff_text(x):
    return x.replace("%2C", ",").replace("%20", " ")


# -----------------------------
# FASTA
# -----------------------------

FASTA = None

for candidate in FASTA_CANDIDATES:
    if os.path.exists(candidate):
        FASTA = candidate
        break

if FASTA is None:
    raise FileNotFoundError("No FASTA found among FASTA_CANDIDATES")

if not os.path.exists(FASTA + ".fai"):
    subprocess.check_call(["samtools", "faidx", FASTA])

fasta_ids = set()

with open(FASTA + ".fai") as handle:
    for line in handle:
        fasta_ids.add(line.split("\t")[0])


# -----------------------------
# CM <-> NC
# -----------------------------

cm_to_nc = {}
nc_to_cm = {}

with open(MAP) as handle:
    for line in handle:
        cm, nc = line.strip().split()[:2]
        cm_to_nc[cm] = nc
        nc_to_cm[nc] = cm


def fasta_chrom_name(gff_chrom):
    if gff_chrom in fasta_ids:
        return gff_chrom

    if gff_chrom in nc_to_cm:
        cm = nc_to_cm[gff_chrom]
        if cm in fasta_ids:
            return cm

    raise ValueError(f"Cannot map {gff_chrom} to FASTA identifiers")


def faidx(chrom, start, end):
    chrom_for_fasta = fasta_chrom_name(chrom)
    region = f"{chrom_for_fasta}:{start}-{end}"

    out = run(["samtools", "faidx", FASTA, region])

    seq = []
    for line in out.splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())

    return "".join(seq).upper()


# -----------------------------
# Read rare CDS variants
# -----------------------------

rare_rows = {}
positions = []

with open(RARE_ANNOT) as handle:
    header = handle.readline().rstrip("\n").split("\t")
    h = {name: i for i, name in enumerate(header)}

    required = ["SNP", "CHROM", "POS", "regulatory_class"]
    for col in required:
        if col not in h:
            raise ValueError(f"Missing column in rare annotation table: {col}")

    for line in handle:
        parts = line.rstrip("\n").split("\t")

        if parts[h["regulatory_class"]] != "CDS":
            continue

        snp = parts[h["SNP"]]
        chrom = parts[h["CHROM"]]
        pos = int(parts[h["POS"]])

        rare_rows[(chrom, pos)] = {
            "SNP": snp,
            "CHROM": chrom,
            "POS": pos,
            "rare_class": parts[h.get("rare_class", -1)] if "rare_class" in h else "NA",
            "AF_TT": parts[h.get("AF_TT", -1)] if "AF_TT" in h else "NA",
            "AF_TX": parts[h.get("AF_TX", -1)] if "AF_TX" in h else "NA",
            "abs_delta_AF": parts[h.get("abs_delta_AF", -1)] if "abs_delta_AF" in h else "NA",
            "nearest_gene_TSS": parts[h.get("nearest_gene_TSS", -1)] if "nearest_gene_TSS" in h else "NA",
            "nearest_gene_description_TSS": parts[h.get("nearest_gene_description_TSS", -1)] if "nearest_gene_description_TSS" in h else "NA",
            "rare_variant_priority_reason": parts[h.get("rare_variant_priority_reason", -1)] if "rare_variant_priority_reason" in h else "NA",
        }

        positions.append((chrom, pos))

if len(positions) == 0:
    with open(OUT, "w") as out:
        out.write("SNP\tCHROM\tPOS\tREF\tALT\tPredicted_effect\n")
    with open(SUMMARY, "w") as out:
        out.write("Predicted_effect\tN\n")
    print("No CDS rare variants found.")
    raise SystemExit(0)


# -----------------------------
# Extract REF / ALT from raw VCF
# -----------------------------

with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
    pos_file = tmp.name
    for chrom, pos in positions:
        tmp.write(f"{chrom}\t{pos}\n")

vcf_cmd = [
    "bcftools",
    "query",
    "-R",
    pos_file,
    "-f",
    "%CHROM\t%POS\t%REF\t%ALT\n",
    RAW_VCF,
]

variant_alleles = {}

for line in run(vcf_cmd).splitlines():
    chrom, pos, ref, alt = line.split("\t")
    pos = int(pos)

    if len(ref) != 1 or len(alt) != 1:
        continue

    variant_alleles[(chrom, pos)] = (ref.upper(), alt.upper())

os.remove(pos_file)


# -----------------------------
# Identify overlapping CDS transcripts
# -----------------------------

target_by_nc = defaultdict(list)

for chrom, pos in positions:
    nc = cm_to_nc.get(chrom, chrom)
    target_by_nc[nc].append((chrom, pos))

parents_to_keep = set()
overlap_records = defaultdict(list)

with open(GFF) as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attr = parts

        if feature != "CDS":
            continue

        if chrom not in target_by_nc:
            continue

        start = int(start)
        end = int(end)

        for cm_chrom, pos in target_by_nc[chrom]:
            if start <= pos <= end:
                parent = parse_attr(attr, "Parent")
                gene = parse_attr(attr, "gene")
                product = decode_gff_text(parse_attr(attr, "product"))
                protein = parse_attr(attr, "protein_id")

                parents_to_keep.add(parent)

                overlap_records[(cm_chrom, pos)].append({
                    "Parent": parent,
                    "Gene_ID": gene,
                    "Product": product,
                    "Protein_ID": protein,
                })


# -----------------------------
# Collect all CDS intervals for retained transcripts
# -----------------------------

cds_by_tx = defaultdict(list)
gene_by_tx = {}
product_by_tx = {}
protein_by_tx = {}
strand_by_tx = {}

with open(GFF) as handle:
    for line in handle:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")

        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attr = parts

        if feature != "CDS":
            continue

        parent = parse_attr(attr, "Parent")

        if parent not in parents_to_keep:
            continue

        gene = parse_attr(attr, "gene")
        product = decode_gff_text(parse_attr(attr, "product"))
        protein = parse_attr(attr, "protein_id")

        cds_by_tx[parent].append((chrom, int(start), int(end), strand, phase))
        gene_by_tx[parent] = gene
        product_by_tx[parent] = product
        protein_by_tx[parent] = protein
        strand_by_tx[parent] = strand


# -----------------------------
# Predict effects
# -----------------------------

rows = []

for tx, intervals in cds_by_tx.items():

    strand = strand_by_tx[tx]

    if strand == "+":
        intervals_ordered = sorted(intervals, key=lambda x: x[1])
    else:
        intervals_ordered = sorted(intervals, key=lambda x: x[2], reverse=True)

    cds_seq_parts = []
    coord_map = []
    offset = 0

    for chrom, start, end, interval_strand, phase in intervals_ordered:
        seq = faidx(chrom, start, end)

        if interval_strand == "-":
            seq = revcomp(seq)

        cds_seq_parts.append(seq)

        length = end - start + 1
        coord_map.append((chrom, start, end, offset, offset + length - 1))
        offset += length

    cds_seq = "".join(cds_seq_parts)

    transcript_chrom = intervals_ordered[0][0]
    cm_chrom = nc_to_cm.get(transcript_chrom, transcript_chrom)

    for chrom, pos in positions:

        if chrom != cm_chrom:
            continue

        if (chrom, pos) not in variant_alleles:
            continue

        ref, alt = variant_alleles[(chrom, pos)]

        snp_offset = None

        for gff_chrom, start, end, off_start, off_end in coord_map:
            if start <= pos <= end:
                if strand == "+":
                    snp_offset = off_start + (pos - start)
                else:
                    snp_offset = off_start + (end - pos)
                break

        if snp_offset is None:
            continue

        ref_coding = ref
        alt_coding = alt

        if strand == "-":
            ref_coding = revcomp(ref_coding)
            alt_coding = revcomp(alt_coding)

        observed_base = cds_seq[snp_offset]
        ref_match = observed_base == ref_coding

        codon_start = (snp_offset // 3) * 3
        codon_ref = cds_seq[codon_start:codon_start + 3]

        if len(codon_ref) != 3:
            continue

        pos_in_codon = snp_offset % 3

        codon_alt_list = list(codon_ref)
        codon_alt_list[pos_in_codon] = alt_coding
        codon_alt = "".join(codon_alt_list)

        aa_ref = CODON_TABLE.get(codon_ref, "X")
        aa_alt = CODON_TABLE.get(codon_alt, "X")

        if aa_ref == aa_alt:
            effect = "synonymous"
        elif aa_alt == "*":
            effect = "stop_gained"
        elif aa_ref == "*":
            effect = "stop_lost"
        else:
            effect = "missense"

        aa_pos = codon_start // 3 + 1

        rr = rare_rows[(chrom, pos)]

        rows.append([
            rr["SNP"],
            chrom,
            pos,
            ref,
            alt,
            rr["rare_class"],
            rr["AF_TT"],
            rr["AF_TX"],
            rr["abs_delta_AF"],
            gene_by_tx[tx],
            product_by_tx[tx],
            tx,
            protein_by_tx[tx],
            strand,
            snp_offset + 1,
            aa_pos,
            codon_ref,
            aa_ref,
            codon_alt,
            aa_alt,
            effect,
            observed_base,
            ref_coding,
            "TRUE" if ref_match else "FALSE",
            rr["nearest_gene_TSS"],
            rr["nearest_gene_description_TSS"],
            rr["rare_variant_priority_reason"],
        ])

header = [
    "SNP",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "rare_class",
    "AF_TT",
    "AF_TX",
    "abs_delta_AF",
    "Gene_ID",
    "Product",
    "Transcript",
    "Protein_ID",
    "Strand",
    "CDS_position_1based",
    "AA_position",
    "Codon_REF",
    "AA_REF",
    "Codon_ALT",
    "AA_ALT",
    "Predicted_effect",
    "Observed_CDS_base",
    "Expected_REF_coding_base",
    "REF_match",
    "nearest_gene_TSS",
    "nearest_gene_description_TSS",
    "rare_variant_priority_reason",
]

with open(OUT, "w") as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")


# -----------------------------
# Summary
# -----------------------------

summary_counts = defaultdict(int)
unique_variant_effect = {}

for row in rows:
    snp = row[0]
    effect = row[20]
    summary_counts[effect] += 1
    unique_variant_effect.setdefault(snp, set()).add(effect)

unique_effect_counts = defaultdict(int)

for snp, effects in unique_variant_effect.items():
    effect_label = ";".join(sorted(effects))
    unique_effect_counts[effect_label] += 1

with open(SUMMARY, "w") as out:
    out.write("summary_type\tPredicted_effect\tN\n")

    for effect, n in sorted(summary_counts.items()):
        out.write(f"transcript_rows\t{effect}\t{n}\n")

    for effect, n in sorted(unique_effect_counts.items()):
        out.write(f"unique_variants\t{effect}\t{n}\n")

print("Done")
print("FASTA:", FASTA)
print("Rare CDS input variants:", len(positions))
print("VCF alleles found:", len(variant_alleles))
print("Overlapping CDS transcripts:", len(parents_to_keep))
print("Predicted rows:", len(rows))
print("Output:", OUT)
print("Summary:", SUMMARY)
