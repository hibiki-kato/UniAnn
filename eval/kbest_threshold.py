#!/usr/bin/env python3
"""Precision / sensitivity of a score threshold on the chrX local k-best paths.

Prediction unit: one distinct proposed CDS structure in one search interval.
Correct: the structure is present in the pCDS reference. Reference structures
are tuples of genomic coordinates and are unique genome-wide, so no locus
restriction is applied; restricting to the locus mapped to the interval would
discard the second gene of a multi-gene path.

Score: the difference from the UniAnn global best path on the same interval, as
reported by the decoder. Raw path scores are not comparable across loci because
the intervals differ in length. A positive difference means the local decoder
outscored the production Viterbi.
"""

import argparse
import csv
import json
from pathlib import Path

from kbest_plotdata import build_predictions

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"           # inputs; gitignored by the repository
OUT = DATA / "eval"            # outputs; gitignored with the rest of data/

CUTS = [500, 0, -100, -250, -500, -750, -1000, -1500, -2000, -3000, -4000,
        -5000, -7500, -10000, None]


def sweep(predictions, total_targets):
    predictions = sorted(predictions, key=lambda r: -r["delta"])
    correct = [r["correct"] for r in predictions]
    deltas = [r["delta"] for r in predictions]
    seen, running = set(), []
    for r in predictions:
        if r["correct"]:
            seen.add(r["structure"])
        running.append(len(seen))

    rows = []
    for cut in CUTS:
        n = len(deltas) if cut is None else sum(1 for d in deltas if d >= cut)
        if n == 0:
            continue
        precision = 100 * sum(correct[:n]) / n
        sensitivity = 100 * running[n - 1] / total_targets
        rows.append({
            "delta_threshold": "none" if cut is None else cut,
            "selected_predictions": n,
            "correct_predictions": sum(correct[:n]),
            "precision_pct": round(precision, 2),
            "structures_recovered": running[n - 1],
            "sensitivity_pct": round(sensitivity, 2),
            "f1_pct": round(2 * precision * sensitivity / (precision + sensitivity), 2)
            if precision + sensitivity else 0.0,
        })

    best = None
    for n in range(1, len(predictions) + 1):
        precision = sum(correct[:n]) / n
        recall = running[n - 1] / total_targets
        if precision + recall == 0:
            continue
        f1 = 2 * precision * recall / (precision + recall)
        if best is None or f1 > best["f1"]:
            best = {"f1": f1, "n": n, "delta": deltas[n - 1],
                    "precision": precision, "recall": recall}
    peak = {
        "delta_threshold": round(best["delta"], 1),
        "selected_predictions": best["n"],
        "precision_pct": round(100 * best["precision"], 2),
        "sensitivity_pct": round(100 * best["recall"], 2),
        "f1_pct": round(100 * best["f1"], 2),
    } if best else None
    return rows, peak, running[-1] if running else 0


def write(rows, path):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def show(title, rows):
    print(title)
    header = list(rows[0])
    print("\t".join(header))
    for r in rows:
        print("\t".join(str(r[h]) for h in header))
    print()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", type=Path, default=DATA / "NC_004354.4.p.CDS.uniq.gtf")
    p.add_argument("--kbest", type=Path, default=DATA / "chrX_local_k_best.k20.gff")
    p.add_argument("--report", type=Path, default=DATA / "chrX_local_k_best.k20.tsv")
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    predictions, truth, total = build_predictions(args)
    all_rows, all_peak, all_recoverable = sweep(predictions, total)

    rank1 = {r["interval"]: r["structure"] for r in predictions if r["rank"] == 1}
    alternatives = [r for r in predictions
                    if rank1.get(r["interval"]) != r["structure"]]
    alt_targets = sum(1 for interval, structures in truth.items()
                      for s in structures if rank1.get(interval) != s)
    alt_rows, alt_peak, alt_recoverable = sweep(alternatives, alt_targets)

    write(all_rows, args.out / "threshold_all.tsv")
    write(alt_rows, args.out / "threshold_alternatives.tsv")

    summary = {
        "score": "difference from the UniAnn global best path on the same interval",
        "all_structures": {
            "predictions": len(predictions),
            "targets": total,
            "recoverable_at_k20": all_recoverable,
            "f1_maximum": all_peak,
        },
        "alternatives_only": {
            "predictions": len(alternatives),
            "targets": alt_targets,
            "recoverable_at_k20": alt_recoverable,
            "f1_maximum": alt_peak,
        },
    }
    (args.out / "threshold_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print()
    show("ALL REFERENCE STRUCTURES", all_rows)
    show("ALTERNATIVES ONLY (rank-1 structure excluded)", alt_rows)


if __name__ == "__main__":
    main()
