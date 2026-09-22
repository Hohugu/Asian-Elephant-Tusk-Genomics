#!/usr/bin/env python3
import argparse
import gzip
import sys
import math

def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]

def gt_to_alt(gt):
    gt = gt.split(":")[0]
    gt = gt.replace("|", "/")
    if gt in (".", "./.", ".|."):
        return None
    parts = gt.split("/")
    if len(parts) != 2:
        return None
    if "." in parts:
        return None
    try:
        a = [int(x) for x in parts]
    except ValueError:
        return None
    if any(x not in (0, 1) for x in a):
        return None
    alt = sum(a)
    het = 1 if alt == 1 else 0
    return alt, het

def calc_group(fields, indices):
    n_called = 0
    alt_sum = 0
    n_het = 0

    for i in indices:
        z = gt_to_alt(fields[4 + i])
        if z is None:
            continue
        alt, het = z
        n_called += 1
        alt_sum += alt
        n_het += het

    if n_called == 0:
        return n_called, alt_sum, float("nan"), float("nan"), float("nan")

    p_alt = alt_sum / (2.0 * n_called)
    hobs = n_het / float(n_called)
    hexp = 2.0 * p_alt * (1.0 - p_alt)

    return n_called, alt_sum, p_alt, hobs, hexp

parser = argparse.ArgumentParser()
parser.add_argument("--sample-order", required=True)
parser.add_argument("--tt", required=True)
parser.add_argument("--tx", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

samples = read_list(args.sample_order)
tt = read_list(args.tt)
tx = read_list(args.tx)

sample_to_idx = {s:i for i,s in enumerate(samples)}

missing_tt = [s for s in tt if s not in sample_to_idx]
missing_tx = [s for s in tx if s not in sample_to_idx]

if missing_tt or missing_tx:
    sys.stderr.write("Missing TT samples: %s\n" % missing_tt[:10])
    sys.stderr.write("Missing TX samples: %s\n" % missing_tx[:10])
    sys.exit(1)

tt_idx = [sample_to_idx[s] for s in tt]
tx_idx = [sample_to_idx[s] for s in tx]

opener = gzip.open if args.out.endswith(".gz") else open

with opener(args.out, "wt") as out:
    out.write("\t".join([
        "CHROM","POS","REF","ALT",
        "TT_n_called","TT_alt_sum","TT_p_alt","TT_hobs","TT_hexp",
        "TX_n_called","TX_alt_sum","TX_p_alt","TX_hobs","TX_hexp",
        "DXY","HetDev_TX","Campbell_Het_candidate"
    ]) + "\n")

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")

        chrom, pos, ref, alt = f[0], f[1], f[2], f[3]

        TT_n, TT_alt, TT_p, TT_hobs, TT_hexp = calc_group(f, tt_idx)
        TX_n, TX_alt, TX_p, TX_hobs, TX_hexp = calc_group(f, tx_idx)

        if math.isnan(TT_p) or math.isnan(TX_p):
            dxy = float("nan")
            hetdev = float("nan")
            camp = "FALSE"
        else:
            dxy = TT_p * (1.0 - TX_p) + TX_p * (1.0 - TT_p)
            hetdev = TX_hobs - TX_hexp
            camp = "TRUE" if (TX_p >= 0.15 and TX_p <= 0.85 and TT_alt == 0) else "FALSE"

        out.write("\t".join(map(str, [
            chrom, pos, ref, alt,
            TT_n, TT_alt, TT_p, TT_hobs, TT_hexp,
            TX_n, TX_alt, TX_p, TX_hobs, TX_hexp,
            dxy, hetdev, camp
        ])) + "\n")
