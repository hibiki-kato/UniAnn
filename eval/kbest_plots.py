#!/usr/bin/env python3
"""Figures for the local k-best score analysis.

Reads the .dat files written by kbest_plotdata.py, so the figures, any LaTeX
plots built from the same files, and the threshold tables all come from one
computation. Requires matplotlib.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "eval"

CORRECT = "#2f6db5"
OTHER = "#c4453c"


def save(fig, args, name):
    fig.savefig(args.out / name, dpi=200)
    if args.also_write:
        args.also_write.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.also_write / name, dpi=200)


def read(path):
    lines = open(path).read().split("\n")
    header = lines[0].split()
    rows = [dict(zip(header, map(float, l.split()))) for l in lines[1:] if l.strip()]
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--out", type=Path, default=DATA)
    p.add_argument("--also-write", type=Path, default=None,
                   help="optional second directory to receive the PNGs")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # --- delta distribution ------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, name, title in (
        (axes[0], "chrX_delta_hist_all.dat", "All 1,783 reference structures"),
        (axes[1], "chrX_delta_hist_alt.dat", "Alternatives only"),
    ):
        rows = read(args.data / name)
        x = [r["centre"] for r in rows]
        ax.bar(x, [r["correct_frac"] for r in rows], width=430,
               color=CORRECT, label="exact match")
        ax.bar(x, [r["other_frac"] for r in rows], width=430,
               color=OTHER, alpha=0.6, label="other")
        ax.set_xlabel(r"$\Delta$ from the UniAnn global path")
        ax.set_ylabel("% of class")
        ax.set_title(title, fontsize=10)
        ax.set_xlim(-15250, 1250)
        ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    save(fig, args, "delta_distribution.png")
    plt.close(fig)

    # --- sensitivity / precision -------------------------------------------
    uniann = read(args.data / "chrX_snpr_uniann.dat")[0]
    rank1 = read(args.data / "chrX_snpr_rank1.dat")[0]
    CEILING = 100 * 1070 / 1783   # one isoform per gene cannot exceed this

    def iso_f1(ax, levels, xmax, ymax):
        """Faint iso-F1 contours: P = f R / (2R - f)."""
        for f in levels:
            xs, ys = [], []
            for i in range(1, 2001):
                r = xmax * i / 2000
                if 2 * r - f <= 1e-9:
                    continue
                pr = f * r / (2 * r - f)
                if pr > ymax:
                    continue
                xs.append(r)
                ys.append(pr)
            if not xs:
                continue
            ax.plot(xs, ys, color="0.72", lw=0.7, zorder=0)
            ax.annotate(f"F1={f:g}", (xs[-1], ys[-1]), color="0.55", fontsize=7,
                        textcoords="offset points", xytext=(-24, 3), zorder=0)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, name, colour, title in (
        (axes[0], "chrX_snpr_all.dat", CORRECT, "All 1,783 reference structures"),
        (axes[1], "chrX_snpr_alt.dat", OTHER, "Alternatives only"),
    ):
        rows = read(args.data / name)
        xmax = max([r["sensitivity"] for r in rows] +
                   ([CEILING] if name.endswith("all.dat") else [])) * 1.12
        ymax = max([r["precision"] for r in rows] +
                   ([uniann["precision"], rank1["precision"]]
                    if name.endswith("all.dat") else [])) * 1.15
        iso_f1(ax, (10, 20, 30, 40, 50) if name.endswith("all.dat") else (5, 10, 15),
               xmax, ymax)
        ax.plot([r["sensitivity"] for r in rows], [r["precision"] for r in rows],
                color=colour, lw=1.8, label="local $K$=20, score threshold")
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, ymax)
        if name.endswith("all.dat"):
            ax.axvline(CEILING, color="0.55", ls="--", lw=0.9, zorder=1)
            ax.annotate(f"one isoform per gene\ncannot exceed {CEILING:.1f}%",
                        (CEILING, ymax * 0.97), fontsize=7.5, color="0.4",
                        ha="right", va="top",
                        textcoords="offset points", xytext=(-5, 0))
            ax.plot(rank1["sensitivity"], rank1["precision"], "D",
                    color="#3f8f5b", markersize=6, markeredgecolor="black",
                    markeredgewidth=0.5, zorder=5, label="local rank-1 only")
            ax.plot(uniann["sensitivity"], uniann["precision"], "*",
                    color="#e08a1e", markersize=15, markeredgecolor="black",
                    markeredgewidth=0.6, zorder=6,
                    label="UniAnn annotation")
            ax.annotate(f"UniAnn {uniann['precision']:.1f}% / {uniann['sensitivity']:.1f}%\n"
                        f"rank-1 {rank1['precision']:.1f}% / {rank1['sensitivity']:.1f}%",
                        (uniann["sensitivity"], uniann["precision"]),
                        textcoords="offset points", xytext=(-30, -34), fontsize=7.5,
                        ha="right")
            ax.legend(frameon=False, fontsize=8, loc="lower left")
        else:
            ax.annotate("UniAnn contributes no\nalternative structures",
                        (0.97, 0.95), xycoords="axes fraction",
                        ha="right", va="top", fontsize=8, color="0.35")
            ax.legend(frameon=False, fontsize=8, loc="lower left")
        marks = read(args.data /
                     ("chrX_snpr_marks_all.dat" if name.endswith("all.dat")
                      else "chrX_snpr_marks_alt.dat"))
        for m in marks:
            ax.plot(m["sensitivity"], m["precision"], "o",
                    color="black", markersize=3.5, zorder=4)
            ax.annotate(f"{m['label']:g}", (m["sensitivity"], m["precision"]),
                        textcoords="offset points", xytext=(5, 4), fontsize=8)
        ax.set_xlabel("sensitivity (%)")
        ax.set_ylabel("precision (%)")
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    save(fig, args, "sensitivity_precision.png")
    plt.close(fig)
    print("wrote", args.out, *(( "and", args.also_write) if args.also_write else ()))


if __name__ == "__main__":
    main()
