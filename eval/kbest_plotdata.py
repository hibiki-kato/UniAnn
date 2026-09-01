#!/usr/bin/env python3
"""Emit plain-text plot data for the local k-best score analysis.

Delta is measured against the UniAnn global best path on the same interval,
which is the quantity reported as score_delta by the decoder. Two views are
written: all reference structures, and alternatives only (the rank-1 structure
of each interval is removed from both the candidates and the targets).
"""

import argparse
import collections
import csv
from pathlib import Path

import kbest_recovery as base

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"           # inputs; gitignored by the repository
OUT = DATA / "eval"            # outputs; gitignored with the rest of data/

BIN_WIDTH = 500
BIN_LO = -15000
BIN_HI = 1000


def build_predictions(args):
    loci, tx_locus, ref = base.read_reference(args.reference)
    kbest_loci, paths = base.read_kbest(args.kbest)
    mapping, _ = base.map_loci_to_genes(loci, kbest_loci)

    # Reference CDS structures are tuples of genomic coordinates and are unique
    # genome-wide, so correctness is membership in the reference set. Restricting
    # the comparison to the locus mapped to the interval would discard the second
    # gene of a multi-gene path, which belongs to a neighbouring locus.
    reference_structures = set(ref.values())
    truth = collections.defaultdict(set)
    for t, structure in ref.items():
        interval = mapping.get(tx_locus[t])
        if interval:
            truth[interval].add(structure)

    delta = {}
    for row in csv.DictReader(open(args.report), delimiter="\t"):
        delta[("locus" + row["reference_gene_index"], int(row["rank"]))] = \
            float(row["score_delta"])

    best = {}
    for interval, entries in paths.items():
        for rank, _gene, structure in entries:
            value = delta.get((interval, rank))
            if value is None:
                continue
            key = (interval, structure)
            if key not in best or value > best[key][0]:
                best[key] = (value, rank, structure in reference_structures)
    predictions = [{"interval": i, "structure": s, "delta": v,
                    "rank": r, "correct": c}
                   for (i, s), (v, r, c) in best.items()]
    predictions.sort(key=lambda r: -r["delta"])
    return predictions, truth, len(ref)


def histogram(predictions, path):
    edges = list(range(BIN_LO, BIN_HI + BIN_WIDTH, BIN_WIDTH))
    pos = collections.Counter()
    neg = collections.Counter()
    for r in predictions:
        index = min(max((r["delta"] - BIN_LO) // BIN_WIDTH, 0), len(edges) - 2)
        (pos if r["correct"] else neg)[int(index)] += 1
    n_pos = max(sum(pos.values()), 1)
    n_neg = max(sum(neg.values()), 1)
    with open(path, "w") as fh:
        fh.write("centre correct other correct_frac other_frac\n")
        for i in range(len(edges) - 1):
            centre = edges[i] + BIN_WIDTH / 2
            fh.write(f"{centre:.0f} {pos[i]} {neg[i]} "
                     f"{100*pos[i]/n_pos:.4f} {100*neg[i]/n_neg:.4f}\n")


def curve(predictions, total_targets, path, step=25, thresholds=()):
    """One point per distinct threshold.

    Predictions that share a delta cannot be separated by a threshold, so a
    point is emitted only at the end of each run of equal deltas. Around
    delta = 0 this matters: roughly 1,100 rank-1 paths tie there, and emitting
    a point inside that run would draw the curve through operating points that
    no threshold can select.
    """
    correct = [r["correct"] for r in predictions]
    deltas = [r["delta"] for r in predictions]
    seen, running = set(), []
    for r in predictions:
        if r["correct"]:
            seen.add(r["structure"])
        running.append(len(seen))
    with open(path, "w") as fh:
        fh.write("sensitivity precision delta selected\n")
        hits = 0
        since = 0
        for n in range(1, len(predictions) + 1):
            hits += correct[n - 1]
            since += 1
            last = n == len(predictions)
            tie = not last and deltas[n] == deltas[n - 1]
            # always emit the exact operating point of a reported threshold
            boundary = not tie and any(
                deltas[n - 1] >= t and (last or deltas[n] < t) for t in thresholds)
            if tie or (since < step and not last and not boundary):
                continue
            since = 0
            fh.write(f"{100*running[n-1]/total_targets:.4f} {100*hits/n:.4f} "
                     f"{deltas[n-1]:.1f} {n}\n")


def uniann_point(reference_path, annotation_path, out_path):
    """Sensitivity and precision of the UniAnn annotation itself.

    UniAnn returns one CDS structure per gene, so it is a single operating
    point rather than a curve.
    """
    def load(path):
        cds = collections.defaultdict(list)
        for line in open(path):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = base.parse_attributes(f[8])
            cds[a["Parent"]].append((int(f[3]), int(f[4]), f[7]))
        return {t: tuple(sorted(v)) for t, v in cds.items()}

    reference = set(load(reference_path).values())
    predicted = set(load(annotation_path).values())
    hits = len(predicted & reference)
    with open(out_path, "w") as fh:
        fh.write("sensitivity precision predictions correct\n")
        fh.write(f"{100*hits/len(reference):.4f} {100*hits/len(predicted):.4f} "
                 f"{len(predicted)} {hits}\n")
    return len(predicted), hits, len(reference)


def marks(curve_path, out_path, thresholds):
    """The points of a curve closest to a set of round threshold values."""
    rows = []
    header = None
    for i, line in enumerate(open(curve_path)):
        if i == 0:
            header = line.split()
            continue
        if line.strip():
            rows.append(dict(zip(header, map(float, line.split()))))
    with open(out_path, "w") as fh:
        fh.write("sensitivity precision label\n")
        for value in thresholds:
            # the operating point of a threshold is the last point that still
            # satisfies it, so that ties at the threshold are all included
            kept = [r for r in rows if r["delta"] >= value]
            if not kept:
                continue
            point = kept[-1]
            fh.write(f"{point['sensitivity']:.4f} {point['precision']:.4f} {value:g}\n")


def rank1_point(predictions, total_targets, out_path):
    """The operating point of taking only the rank-1 path of every interval.

    This is the like-for-like comparison with UniAnn: one decoded path per
    search interval. It is not the same as thresholding at delta >= 0, which
    also admits lower-ranked paths from the intervals where the local decoder
    outscores the global one.
    """
    rows = [r for r in predictions if r["rank"] == 1]
    hits = len({r["structure"] for r in rows if r["correct"]})
    correct = sum(r["correct"] for r in rows)
    with open(out_path, "w") as fh:
        fh.write("sensitivity precision predictions correct\n")
        fh.write(f"{100*hits/total_targets:.4f} {100*correct/len(rows):.4f} "
                 f"{len(rows)} {correct}\n")
    return len(rows), correct, hits


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", type=Path, default=DATA / "NC_004354.4.p.CDS.uniq.gtf")
    p.add_argument("--kbest", type=Path, default=DATA / "chrX_local_k_best.k20.gff")
    p.add_argument("--report", type=Path, default=DATA / "chrX_local_k_best.k20.tsv")
    p.add_argument("--annotation", type=Path,
                   default=DATA / "NC_004354.4.fa.uniann.gff",
                   help="the ordinary UniAnn annotation, plotted as one operating point")
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    predictions, truth, total = build_predictions(args)
    marked = (0, -500, -1000, -2000, -5000)
    histogram(predictions, args.out / "chrX_delta_hist_all.dat")
    curve(predictions, total, args.out / "chrX_snpr_all.dat", thresholds=marked)

    rank1 = {r["interval"]: r["structure"] for r in predictions if r["rank"] == 1}
    alternatives = [r for r in predictions
                    if rank1.get(r["interval"]) != r["structure"]]
    alt_targets = sum(1 for interval, structures in truth.items()
                      for s in structures if rank1.get(interval) != s)
    histogram(alternatives, args.out / "chrX_delta_hist_alt.dat")
    curve(alternatives, alt_targets, args.out / "chrX_snpr_alt.dat",
          thresholds=marked)

    marks(args.out / "chrX_snpr_all.dat", args.out / "chrX_snpr_marks_all.dat", marked)
    marks(args.out / "chrX_snpr_alt.dat", args.out / "chrX_snpr_marks_alt.dat", marked)

    n_r1, c_r1, h_r1 = rank1_point(predictions, total,
                                   args.out / "chrX_snpr_rank1.dat")
    print(f"local rank-1 only: {n_r1} structures, {c_r1} correct, "
          f"precision {100*c_r1/n_r1:.2f}%, sensitivity {100*h_r1/total:.2f}%")

    n_pred, n_hit, n_ref = uniann_point(
        args.reference, args.annotation,
        args.out / "chrX_snpr_uniann.dat")
    print(f"UniAnn annotation: {n_pred} structures, {n_hit} correct, "
          f"precision {100*n_hit/n_pred:.2f}%, sensitivity {100*n_hit/n_ref:.2f}%")
    print(f"all: {len(predictions)} predictions, {total} targets")
    print(f"alt: {len(alternatives)} predictions, {alt_targets} targets")
    for name in ("chrX_delta_hist_all.dat", "chrX_delta_hist_alt.dat",
                 "chrX_snpr_all.dat", "chrX_snpr_alt.dat"):
        print(f"  {name}: {sum(1 for _ in open(args.out / name)) - 1} rows")


if __name__ == "__main__":
    main()
