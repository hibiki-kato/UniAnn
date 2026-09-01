# Evaluating local k-best output

Python tools that score the annotation and the `--local-k-best` output of
`uniann` against a reference annotation. They are not part of the build and
`install.sh` does not copy them; run them from this directory.

`kbest_recovery.py`, `kbest_threshold.py` and `kbest_plotdata.py` use the
standard library alone. `kbest_plots.py` needs matplotlib and
`plot_best_path_score.py` needs numpy, pandas and matplotlib.

This branch is the decoder branch plus this directory. Nothing here is imported
by `uniann`, and `install.sh` does not copy it; keeping it off the decoder
branch keeps that branch free of Python dependencies.

## Inputs

- a reference annotation with `locus`, `transcript` and `CDS` records on the
  forward strand, for example the output of `gffread -M` on a curated set
- `--local-k-output` GFF3 and `--local-k-report` TSV from a `uniann` run
- the ordinary UniAnn annotation, used as a single operating point

Paths default to `../data/`, which the repository ignores. Every tool takes
explicit paths, so nothing here assumes a particular experiment.

## Tools

| script | what it does | writes |
|---|---|---|
| `kbest_recovery.py` | exact phase-aware CDS recovery by rank and by subset | `recovery_transcripts.tsv`, `recovery_summary.json` |
| `kbest_threshold.py` | precision and sensitivity of a score threshold | `threshold_all.tsv`, `threshold_alternatives.tsv`, `threshold_summary.json` |
| `kbest_plotdata.py` | plot data for the two views, plus the UniAnn and rank-1 operating points | `*.dat` |
| `kbest_plots.py` | renders the `.dat` files | `delta_distribution.png`, `sensitivity_precision.png` |
| `plot_best_path_score.py` | plots the Viterbi best-path score along the sequence from the `out.err` DP trace | a PNG per requested window |

Output defaults to `../data/eval/`, also ignored. `--out` moves it, and
`kbest_plots.py --also-write DIR` copies the PNGs somewhere else as well.

`kbest_threshold.py` imports `kbest_plotdata.py`, which imports
`kbest_recovery.py`, so the tables, the plot data and the figures come from one
computation and cannot disagree.

## Definitions

Predictions are deduplicated to distinct CDS structures. A prediction is
correct when its structure is present in the reference; reference CDS
structures are tuples of genomic coordinates and are unique genome-wide, so no
locus restriction is applied. Restricting the comparison to the locus mapped to
a search interval would discard the second gene of a multi-gene path.

The threshold is applied to the score difference from the UniAnn global best
path on the same interval. Raw path scores are not comparable across loci
because the intervals differ in length. A positive difference means the local
decoder outscored the production Viterbi.

Predictions that share a value of the difference cannot be separated by a
threshold, so each curve carries one point per distinct value.

## Example

```
uniann.sh -f chrX.fa -s sites.tsv -p psauron.csv -n \
          --local-k-best 20 \
          --local-k-output kbest.gff --local-k-report kbest.tsv

cd eval
python3 kbest_recovery.py  --reference ref.gtf --kbest ../kbest.gff
python3 kbest_threshold.py --reference ref.gtf --kbest ../kbest.gff --report ../kbest.tsv
python3 kbest_plotdata.py  --reference ref.gtf --kbest ../kbest.gff --report ../kbest.tsv \
                           --annotation ../chrX.fa.uniann.gff
python3 kbest_plots.py
```
