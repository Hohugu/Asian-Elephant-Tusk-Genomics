#!/usr/bin/env python3

import argparse
import bisect
import csv
import gzip
from pathlib import Path
from collections import Counter, defaultdict

METRICS = ["FST", "DXY", "muLD", "HetDev_TX"]

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = []

    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e

    merged.append((cs, ce))
    return merged

def load_top_windows(top_windows_file):
    raw = {m: defaultdict(list) for m in METRICS}

    with open_text(top_windows_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            metric = row["metric"]
            if metric not in raw:
                continue

            chrom = row["CHROM"]
            start = int(float(row["window_start"]))
            end = int(float(row["window_end"]))

            if end < start:
                start, end = end, start

            raw[metric][chrom].append((start, end))

    merged = {m: {} for m in METRICS}
    summary = []

    for metric in METRICS:
        for chrom, intervals in raw[metric].items():
            merged_intervals = merge_intervals(intervals)
            starts = [x[0] for x in merged_intervals]
            ends = [x[1] for x in merged_intervals]

            merged[metric][chrom] = (starts, ends)

            summary.append({
                "metric": metric,
                "CHROM": chrom,
                "n_input_windows": len(intervals),
                "n_merged_intervals": len(merged_intervals),
                "total_bp_merged": sum(e - s + 1 for s, e in merged_intervals)
            })

    return merged, summary

def in_intervals(intervals_by_metric, metric, chrom, pos):
    data = intervals_by_metric.get(metric, {}).get(chrom)
    if data is None:
        return False

    starts, ends = data
    i = bisect.bisect_right(starts, pos) - 1

    if i < 0:
        return False

    return ends[i] >= pos

def process_one_threshold(label, top_windows_file, universe_file, outdir):
    print(f"[{label}] Loading top windows: {top_windows_file}", flush=True)

    intervals, interval_summary = load_top_windows(top_windows_file)

    outdir.mkdir(parents=True, exist_ok=True)

    interval_summary_file = outdir / f"Campbell_SNP_interval_summary_{label}.tsv"
    with open(interval_summary_file, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=["metric", "CHROM", "n_input_windows", "n_merged_intervals", "total_bp_merged"]
        )
        writer.writeheader()
        for row in interval_summary:
            writer.writerow(row)

    membership_file = outdir / f"Campbell_SNP_membership_{label}.tsv.gz"
    counts_file = outdir / f"Campbell_SNP_counts_{label}.tsv"
    pattern_file = outdir / f"Campbell_SNP_pattern_counts_{label}.tsv"
    pairwise_file = outdir / f"Campbell_SNP_pairwise_counts_{label}.tsv"

    metric_counts = Counter()
    pattern_counts = Counter()
    n_analyses_counts = Counter()
    pairwise_counts = Counter()

    total_snps = 0
    selected_snps = 0

    print(f"[{label}] Streaming SNPs from: {universe_file}", flush=True)

    with open_text(universe_file) as f, gzip.open(membership_file, "wt", newline="") as out:
        reader = csv.DictReader(f, delimiter="\t")

        fields = [
            "snp_id", "CHROM", "POS", "REF", "ALT",
            "FST", "DXY", "muLD", "HetDev_TX",
            "n_analyses", "pattern",
            "DXY_value", "HetDev_TX_value",
            "TX_p_alt", "TT_alt_sum", "Campbell_Het_candidate"
        ]

        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()

        for row in reader:
            total_snps += 1

            chrom = row["CHROM"]
            pos = int(float(row["POS"]))

            flags = {}
            active = []

            for metric in METRICS:
                hit = in_intervals(intervals, metric, chrom, pos)
                flags[metric] = hit
                if hit:
                    active.append(metric)

            if not active:
                continue

            selected_snps += 1

            for metric in active:
                metric_counts[metric] += 1

            for i in range(len(active)):
                for j in range(i + 1, len(active)):
                    pairwise_counts[(active[i], active[j])] += 1

            n_analyses = len(active)
            n_analyses_counts[n_analyses] += 1

            pattern = "&".join(active)

            pattern_counts[pattern] += 1

            ref = row.get("REF", "")
            alt = row.get("ALT", "")
            snp_id = f"{chrom}:{row['POS']}:{ref}:{alt}"

            writer.writerow({
                "snp_id": snp_id,
                "CHROM": chrom,
                "POS": row["POS"],
                "REF": ref,
                "ALT": alt,
                "FST": str(flags["FST"]),
                "DXY": str(flags["DXY"]),
                "muLD": str(flags["muLD"]),
                "HetDev_TX": str(flags["HetDev_TX"]),
                "n_analyses": n_analyses,
                "pattern": pattern,
                "DXY_value": row.get("DXY", ""),
                "HetDev_TX_value": row.get("HetDev_TX", ""),
                "TX_p_alt": row.get("TX_p_alt", ""),
                "TT_alt_sum": row.get("TT_alt_sum", ""),
                "Campbell_Het_candidate": row.get("Campbell_Het_candidate", "")
            })

            if total_snps % 1000000 == 0:
                print(f"[{label}] processed {total_snps:,} SNPs; retained {selected_snps:,}", flush=True)

    with open(counts_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["threshold", "category", "count"])

        writer.writerow([label, "total_universe_snps", total_snps])
        writer.writerow([label, "SNPs_in_any_analysis", selected_snps])

        for metric in METRICS:
            writer.writerow([label, f"{metric}_SNPs", metric_counts[metric]])

        for n in sorted(n_analyses_counts):
            writer.writerow([label, f"SNPs_in_{n}_analyses", n_analyses_counts[n]])

    with open(pattern_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["threshold", "pattern", "count"])

        for pattern, count in pattern_counts.most_common():
            writer.writerow([label, pattern, count])

    with open(pairwise_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["threshold", "metric1", "metric2", "count"])

        for i in range(len(METRICS)):
            for j in range(i + 1, len(METRICS)):
                m1 = METRICS[i]
                m2 = METRICS[j]
                writer.writerow([label, m1, m2, pairwise_counts[(m1, m2)]])

    print(f"[{label}] Done.", flush=True)
    print(f"[{label}] Output: {membership_file}", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--camp", required=True)
    args = parser.parse_args()

    camp = Path(args.camp)

    universe_file = camp / "metrics/dxy_het/DXY_HetDev_TX.full.tsv.gz"
    outdir = camp / "snp_level_campbell"

    thresholds = [
        (
            "5pct_Campbell",
            camp / "venn_campbell_corrected/top_windows_5pct_Campbell.tsv.gz"
        ),
        (
            "1pct_Campbell_threshold",
            camp / "venn_campbell_corrected/top_windows_1pct_Campbell_threshold.tsv.gz"
        )
    ]

    for label, top_windows_file in thresholds:
        process_one_threshold(label, top_windows_file, universe_file, outdir)

if __name__ == "__main__":
    main()
