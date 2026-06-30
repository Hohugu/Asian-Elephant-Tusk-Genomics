#!/usr/bin/env python3

import os
import re
import subprocess
from collections import defaultdict

BASE = "/scratch/project_2000886/Hoedric/GWAS_2025"

GFF = "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.gff"

FASTA_CANDIDATES = [
    "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fna",
    "/scratch/project_2000886/Elisa/reseq/GCA_024166365.1_mEleMax1_primary_haplotype_genomic.fa",
    "/scratch/project_2000886/Elisa/reseq/GCF_024166365.1_mEleMax1_primary_haplotype_genomic.fa",
]

MAP = f"{BASE}/GWAS/tables_GEMMA_male/CM_to_NC.map"
VCF = f"{BASE}/Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_CDS_SNPs.raw.vcf.gz"
CDS_TABLE = f"{BASE}/Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_CDS_SNPs.tsv"
OUT = f"{BASE}/Genetics_Analysis/Functional_annotation/regulatory/GWAS_Bonferroni_CDS_SNP_effects.tsv"

CODON = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMP)[::-1].upper()


def attr_value(attr, key):
    m = re.search(rf"{key}=([^;]+)", attr)
    return m.group(1) if m else ""


def clean_text(x):
    return x.replace("%2C", ",").replace("%20", " ")


def run(cmd):
    return subprocess.check_output(cmd, text=True)


# ---------- FASTA ----------
FASTA = next((f for f in FASTA_CANDIDATES if os.path.exists(f)), None)
if FASTA is None:
    raise FileNotFoundError("No FASTA found among FASTA_CANDIDATES")

if not os.path.exists(FASTA + ".fai"):
    subprocess.check_call(["samtools", "faidx", FASTA])

with open(FASTA + ".fai") as f:
    fasta_ids = {line.split("\t")[0] for line in f}


# ---------- CM/NC map ----------
cm_to_nc = {}
nc_to_cm = {}
with open(MAP) as f:
    for line in f:
        cm, nc = line.strip().split()[:2]
        cm_to_nc[cm] = nc
        nc_to_cm[nc] = cm


def chrom_for_fasta(chrom):
    if chrom in fasta_ids:
        return chrom
    if chrom in nc_to_cm and nc_to_cm[chrom] in fasta_ids:
        return nc_to_cm[chrom]
    raise ValueError(f"Cannot map chromosome {chrom} to FASTA identifiers")


def faidx(chrom, start, end):
    chrom_fa = chrom_for_fasta(chrom)
    region = f"{chrom_fa}:{start}-{end}"
    txt = run(["samtools", "faidx", FASTA, region])
    return "".join(x.strip() for x in txt.splitlines() if not x.startswith(">")) .upper()


# ---------- CDS SNP table ----------
genes = set()
models_by_site = defaultdict(set)

with open(CDS_TABLE) as f:
    header = f.readline().rstrip("\n").split("\t")
    h = {name: i for i, name in enumerate(header)}

    for line in f:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[h["CHROM"]]
        pos = fields[h["POS"]]
        site = f"{chrom}:{pos}"

        if "GWAS_model" in h:
            models_by_site[site].add(fields[h["GWAS_model"]])

        if "Gene_ID" in h:
            gene = fields[h["Gene_ID"]]
            if gene not in {"", "NA", "<NA>"}:
                genes.add(gene)

if not genes:
    genes = {"LOC126078159", "LOC126079052"}


# ---------- VCF SNPs ----------
snps = []
vcf_cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", VCF]
for line in run(vcf_cmd).splitlines():
    chrom, pos, ref, alt = line.split("\t")
    if len(ref) != 1 or len(alt) != 1:
        continue
    pos = int(pos)
    snps.append({
        "CM": chrom,
        "NC": cm_to_nc.get(chrom, chrom),
        "pos": pos,
        "ref": ref.upper(),
        "alt": alt.upper(),
        "site": f"{chrom}:{pos}",
    })


# ---------- GFF CDS ----------
cds_by_tx = defaultdict(list)
gene_by_tx = {}
product_by_tx = {}
protein_by_tx = {}
strand_by_tx = {}

with open(GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue

        chrom, source, feature, start, end, score, strand, phase, attr = fields
        gene = attr_value(attr, "gene")
        if gene not in genes:
            continue

        tx = attr_value(attr, "Parent")
        product = clean_text(attr_value(attr, "product"))
        protein = attr_value(attr, "protein_id")

        start = int(start)
        end = int(end)
        cds_by_tx[tx].append((chrom, start, end, strand, phase))
        gene_by_tx[tx] = gene
        product_by_tx[tx] = product
        protein_by_tx[tx] = protein
        strand_by_tx[tx] = strand


# ---------- Effect prediction ----------
rows = []

for tx, intervals in cds_by_tx.items():
    tx_strand = strand_by_tx[tx]

    if tx_strand == "+":
        intervals = sorted(intervals, key=lambda x: x[1])
    else:
        intervals = sorted(intervals, key=lambda x: x[2], reverse=True)

    cds_seq_parts = []
    coord_map = []
    offset = 0

    for chrom, start, end, strand, phase in intervals:
        seq = faidx(chrom, start, end)
        if strand == "-":
            seq = revcomp(seq)

        cds_seq_parts.append(seq)
        length = end - start + 1
        coord_map.append((chrom, start, end, offset, offset + length - 1))
        offset += length

    cds_seq = "".join(cds_seq_parts)

    for snp in snps:
        if snp["NC"] != intervals[0][0]:
            continue

        snp_offset = None
        for chrom, start, end, off_start, off_end in coord_map:
            if start <= snp["pos"] <= end:
                if tx_strand == "+":
                    snp_offset = off_start + (snp["pos"] - start)
                else:
                    snp_offset = off_start + (end - snp["pos"])
                break

        if snp_offset is None:
            continue

        ref_coding = snp["ref"]
        alt_coding = snp["alt"]
        if tx_strand == "-":
            ref_coding = revcomp(ref_coding)
            alt_coding = revcomp(alt_coding)

        observed_base = cds_seq[snp_offset]
        ref_match = observed_base == ref_coding

        codon_start = (snp_offset // 3) * 3
        codon_ref = cds_seq[codon_start:codon_start + 3]
        if len(codon_ref) != 3:
            continue

        pos_in_codon = snp_offset % 3
        codon_alt = list(codon_ref)
        codon_alt[pos_in_codon] = alt_coding
        codon_alt = "".join(codon_alt)

        aa_ref = CODON.get(codon_ref, "X")
        aa_alt = CODON.get(codon_alt, "X")

        if aa_ref == aa_alt:
            effect = "synonymous"
        elif aa_alt == "*":
            effect = "stop_gained"
        elif aa_ref == "*":
            effect = "stop_lost"
        else:
            effect = "missense"

        aa_pos = codon_start // 3 + 1

        rows.append([
            snp["CM"],
            snp["pos"],
            snp["ref"],
            snp["alt"],
            ",".join(sorted(models_by_site.get(snp["site"], []))) or "NA",
            gene_by_tx[tx],
            tx,
            protein_by_tx[tx],
            product_by_tx[tx],
            tx_strand,
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
        ])


# ---------- Output ----------
header = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "GWAS_models",
    "Gene_ID",
    "Transcript",
    "Protein_ID",
    "Product",
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
]

with open(OUT, "w") as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")

print("Done")
print("FASTA:", FASTA)
print("Genes:", ",".join(sorted(genes)))
print("SNPs:", len(snps))
print("Transcripts:", len(cds_by_tx))
print("Predicted rows:", len(rows))
print("Output:", OUT)
