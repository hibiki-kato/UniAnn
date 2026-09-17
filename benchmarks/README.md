# Comparing the two best-path decoders

Build the decoder and place prepared inputs in an ignored local directory, for
example `data/dmel/` or `data/hsap/`. The script needs Python 3.8+, GNU
`/usr/bin/time`, and the repository's `gffread`. It uses no Python packages.

For the default Dmel layout:

```text
data/dmel/
  NC_004354.4.fa
  NC_004354.4.p.CDS.uniq.gtf
  out.ps.txt
  out.gt.txt
  out.ag.txt
  out.atg.txt
  out.stop.txt
```

```sh
make -C src
python3 benchmarks/compare_viterbi.py \
  --data data/dmel --out data/dmel/comparison
```

By default both runs use `src/uniann`: one with `--no-dp-dump`, the other with
that flag plus `--metadata-state-viterbi`. `--baseline-binary` and
`--metadata-binary` can select separately built binaries. The baseline must
support `--no-dp-dump`; pristine upstream does not support that flag.

For another dataset specify `--fasta` and `--reference` explicitly. The five
score filenames above remain relative to `--data`. Use a single-record FASTA
and a reference covering exactly the decoded sequence and strand. The script
does not subset the reference or transform strand coordinates. It accepts
GFF3 `Parent` or GTF `transcript_id` attributes on CDS rows; missing transcript
identifiers cause an error.

## Outputs and interpretation

Each variant directory contains its raw GFF, filtered annotation, decoder and
`gffread` stderr logs, GNU time log, metrics and a manifest. The manifest records
the command, binary/tool SHA-256 hashes, input/reference hashes and filter
threshold. Paired summaries reject differences in the inputs, tool or filter.
Runs should use inputs that remain unchanged throughout the comparison.
Hashing, reference evaluation and `gffread` are outside the timed decoder call.

Postprocessing matches the wrapper: `gffread --tlf`, removal of one-/two-exon
predictions whose summed exon length is at most 200 bases, then conversion back
to GFF. Raw decoder predictions contain CDS exons, so this is the CDS-length
filter used by that workflow. `--min-cds` changes the threshold for both runs.

Exact matching groups CDS rows by sequence, strand and transcript identifier,
then represents each transcript as its sorted `(start, end, phase)` segments.
Identical structures from different transcript identifiers are deduplicated.
A prediction matches only when sequence, strand and the complete segment list
agree; this is a phase-aware exact CDS metric, not gffcompare transcript F1.

- SN = 100 × TP / reference structures.
- PR = 100 × TP / predicted structures.
- F1 = 200 × TP / (reference structures + predicted structures).

`summary.json` includes integer counts, resource measurements and structures
added or removed between variants. `comparison.tsv` contains SN/PR/F1 percentages.
Resource units are seconds and KiB, as named in the JSON fields. Empty
prediction/reference denominators produce zero, not an undefined number.

To rerun one variant, use `--variants metadata` or `--variants baseline`. To
regenerate summaries from successful saved runs, pass `--variants` with no
values and the same reference arguments. A failed rerun removes the old
variant metrics so they cannot silently appear as a new successful run.

See [the feature document](../docs/metadata_viterbi.md#measurements) for measured
Dmel accuracy and memory, and full human chr1 metadata memory. These local data
and generated experiment artifacts remain ignored and are not part of the PR.
