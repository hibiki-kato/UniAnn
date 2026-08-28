#!/usr/bin/env python3
"""Plot the Viterbi best-path score of uniann (out.err dp/bt trace) against genome position.

Reproduces viterbi_termination + viterbi_backtrace_scores from src/uniann.cpp:
  best_state = argmax_s dp[L-1][s]
  path_scores[i] = dp[i][cur]; cur = bt[i][cur]   (walking i = L-1 .. 0)

Usage:
  plot_best_path_score.py DP_TXT BT_TXT OUTDIR [--start S --end E]
where DP_TXT / BT_TXT hold only the 7 numeric state columns per position, e.g.
  awk '$2=="dp"{print $3,$4,$5,$6,$7,$8,$9}' out.err > dp.txt
  awk '$2=="bt"{print $3,$4,$5,$6,$7,$8,$9}' out.err > bt.txt
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NUM_STATES = 7


def load_matrix(path, dtype):
    df = pd.read_csv(path, sep=" ", header=None, dtype=np.int64,
                     names=list(range(NUM_STATES)))
    return df.to_numpy().astype(dtype)


def backtrace(dp, bt):
    L = dp.shape[0]
    states = np.empty(L, dtype=np.int8)
    scores = np.empty(L, dtype=np.int64)
    cur = int(np.argmax(dp[L - 1]))
    for i in range(L - 1, -1, -1):
        states[i] = cur
        scores[i] = dp[i, cur]
        cur = int(bt[i, cur])
        if cur < 0:  # position 0 has bt == -1
            cur = 0
    return states, scores


def load_genes(gff, lo, hi):
    """Per-transcript CDS structure, 0-based half-open dp indices.

    Returns [{id, strand, start, end, exons}] where start/end are the coding
    start (ATG) and stop end of the transcript, i.e. the span the score
    difference scores[end-1] - scores[start-1] refers to. Introns are the gaps
    between consecutive `exons`.
    """
    tx = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            s, e = int(f[3]) - 1, int(f[4])  # GFF is 1-based inclusive
            parent = ""
            for kv in f[8].split(";"):
                if kv.startswith("Parent="):
                    parent = kv[7:]
                    break
            rec = tx.setdefault(parent, {"id": parent, "strand": f[6], "exons": []})
            rec["exons"].append((s, e))
    genes = []
    for rec in tx.values():
        rec["exons"].sort()
        rec["start"] = rec["exons"][0][0]
        rec["end"] = rec["exons"][-1][1]
        if rec["end"] > lo and rec["start"] < hi:
            genes.append(rec)
    genes.sort(key=lambda g: g["start"])
    return genes


def draw_genes(ax, genes):
    """Light band = whole gene (start codon -> stop codon), dark blocks = coding exons."""
    for k, g in enumerate(genes):
        ax.axvspan(g["start"], g["end"], color="#9467bd", alpha=0.10, lw=0, zorder=0,
                   label="gene (start->stop)" if k == 0 else None)
        for j, (s, e) in enumerate(g["exons"]):
            ax.axvspan(s, e, color="#d62728", alpha=0.22, lw=0, zorder=1,
                       label="CDS exon" if k == 0 and j == 0 else None)
        ax.axvline(g["start"], color="#2ca02c", lw=0.8, zorder=2,
                   label="coding start" if k == 0 else None)
        ax.axvline(g["end"], color="#8c564b", lw=0.8, zorder=2,
                   label="stop end" if k == 0 else None)


def envelope(x, y, n_bins):
    """min/max per bin so a decimated line keeps the shape of the full trace."""
    L = len(y)
    if L <= n_bins * 2:
        return x, y
    per = L // n_bins
    trimmed = y[: per * n_bins].reshape(n_bins, per)
    xt = x[: per * n_bins].reshape(n_bins, per)
    lo, hi = trimmed.min(axis=1), trimmed.max(axis=1)
    xi = xt[:, 0]
    # interleave min/max so plt.plot draws a filled-looking envelope
    xs = np.repeat(xi, 2)
    ys = np.empty(2 * n_bins, dtype=y.dtype)
    ys[0::2], ys[1::2] = lo, hi
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dp_txt", help="or, with --from-npy, the saved best_path_scores.npy")
    ap.add_argument("bt_txt", nargs="?")
    ap.add_argument("outdir")
    ap.add_argument("--from-npy", action="store_true",
                    help="skip parsing/backtrace, reuse best_path_scores.npy")
    ap.add_argument("--cds-gff", help="GFF whose CDS features are shaded in the background")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None)
    ap.add_argument("--bins", type=int, default=20000)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if args.from_npy:
        scores = np.load(args.dp_txt)
    else:
        dp = load_matrix(args.dp_txt, np.int64)
        bt = load_matrix(args.bt_txt, np.int8)
        assert dp.shape == bt.shape, (dp.shape, bt.shape)

        states, scores = backtrace(dp, bt)
        del dp, bt
        np.save(os.path.join(args.outdir, "best_path_scores.npy"), scores)
        np.save(os.path.join(args.outdir, "best_path_states.npy"), states)

    end = args.end if args.end is not None else len(scores)
    pos = np.arange(args.start, end, dtype=np.int64)
    sub = scores[args.start:end]

    xs, ys = envelope(pos, sub, args.bins)
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(xs, ys, lw=0.8, color="#1f77b4", zorder=3)
    if args.cds_gff:
        genes = load_genes(args.cds_gff, args.start, end)
        print(f"{len(genes)} genes in window")
        draw_genes(ax, genes)
        if genes:
            ax.legend(loc="lower right", fontsize=8)
    ax.set_xlabel("genomic position (base)")
    ax.set_ylabel("Viterbi best-path score")
    ax.set_title(f"uniann best-path score, {args.start}-{end} "
                 f"({len(sub)} bases, {args.bins} bins min/max)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(args.outdir, f"best_path_score_{args.start}_{end}.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)
    print(f"L={len(scores)} final={scores[-1]} min={sub.min()} max={sub.max()}")


if __name__ == "__main__":
    main()
