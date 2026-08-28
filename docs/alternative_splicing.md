# Experimental anchored alternative-splicing mode

UniAnn can optionally search for a locally optimal complete transcript around
each donor and acceptor selected by its normal global Viterbi path. The normal
prediction is still written to stdout exactly as before. Alternative mode writes
an additional GFF3 file and, optionally, a diagnostic TSV.

```bash
src/uniann sequence.fa emissions.txt gt.txt ag.txt atg.txt stop.txt \
  --alternative-splicing \
  --alt-output prediction.alternative_splicing.gff \
  --alt-report prediction.alternative_splicing.tsv \
  --alt-max-distance 1000000
```

The `scripts/uniann.sh` wrapper accepts the same long options. When its
`--alt-output` is omitted it writes
`<FASTA>.uniann.alternative_splicing.gff`; a direct executable invocation
defaults to `<FASTA>.alternative_splicing.gff`.

Options:

- `--alternative-splicing`: enable the experimental search.
- `--alt-output FILE`: set the additional primary GFF3 output.
- `--alt-report FILE`: write the secondary diagnostic TSV.
- `--alt-max-distance BP`: bound each direction; default 1,000,000.
- `--alt-min-score-delta SCORE`: omit complete unique alternatives below the
  threshold from primary GFF3. The default is no threshold.
- `--alt-include-incomplete`: retain incomplete searches in the TSV only.
- `--alt-debug`: emit per-anchor interval, endpoint, score, failure, and timing
  diagnostics to stderr.

## Anchor and coordinate definitions

Every `E* -> I*` transition on the global path is a donor anchor and every
`I* -> E*` transition is an acceptor anchor. The anchored local path must contain
the same transition states at the same DP position exactly once.

All internal/report `_0based` positions are zero-based. An anchor genomic
position is the first base of its GT or AG dinucleotide. Its `_1based` value is
that position plus one. Intervals are zero-based half-open. GFF3 output is
one-based inclusive.

## Search and scoring

The upstream search is a constrained forward DP, chosen because UniAnn's exon,
intron, and intergenic lengths plus `exon_from` metadata are not recoverable from
a state label during a simple reverse recurrence. A virtual intergenic source
represents all valid starts in one bounded pass. The anchored transition is
forced, including its full metadata. The same pass continues downstream and
compares every valid exon-to-intergenic transcript endpoint in the interval.

Each local cell has reference-splice and splice-diverged hypotheses. The latter
is needed because conditioning an additive Viterbi search only on a transition
already present in the global optimum would otherwise reconstruct the global
path by optimal substructure. Divergence changes only splice transitions; no
score perturbation is used. This returns at most one best divergent candidate
per anchor and is not a k-best algorithm.

The DP score is replayed through the same authoritative transition evaluator.
The candidate is rejected unless the scale-aware score difference is at most
`1e-8`. Validation also checks contiguity, transition/emission legality, motif
coordinates, start/stop codons, frame, segment lengths, the exact anchor,
transcript boundaries, and unambiguous annotation geometry.

`candidate_total_score` is the additive transcript score from its start
transition through its terminating transition. `reference_interval_score`
replays the global path across the same genomic interval, including intergenic
bases after an earlier reference endpoint. `score_delta` is candidate minus
that comparable-interval score. This definition prevents a missing suffix from
being treated as free, but comparisons can still depend on the candidate's
chosen interval.

## Output and determinism

The additional output is valid GFF3 with `gene`, `transcript`, `exon`, `CDS`, and
`intron` features. The original transcript is `.t1`; unique alternatives are
ordered by score then canonical coordinates and named `.t2`, `.t3`, and so on.
All transcripts from an anchored locus share `UniAnnGeneNNNNNN`. CDS phase is
calculated once by the shared transcript writer.

Candidates are primarily deduplicated by sequence, strand, reference locus, and
ordered intron chain. The highest-scoring complete representative is retained.
A structure whose intron chain equals its reference is diagnostic only and is
not emitted again. Classifications are deliberately conservative:
`alternative_donor`, `alternative_acceptor`, `exon_skipping`,
`alternative_exon`, `intron_retention`, or `complex`.

Equal local scores preserve loop order; retained candidates then use descending
score, lexicographic intron chain, genomic start, genomic end, and anchor
coordinate. Ordered containers are used for scientific output.

## Complexity

Let `L` be sequence length, `S` the seven HMM states, `A` the number of anchors,
`R` the maximum bounded interval length, and `C` the retained candidate count.

- Global time: `O(L S^2)`; the preserved legacy loop visits all state pairs.
- Local time: `O(A R S^2)` with a constant factor of four for anchor-seen and
  splice-diverged hypotheses.
- Global memory/backpointers: `O(L S)`.
- One local DP/backpointer buffer: `O(R S)` with the four-hypothesis factor.
- Stored candidate paths/transcripts: `O(C R)`; raw pre-dedup diagnostics can
  additionally retain `O(A R)` path data.

Searches are sequential, so local DP memory is not multiplied by `A`. Immutable
sequence and score arrays are referenced, never copied per anchor.

## Known limitations

- A candidate must share at least one exact donor or acceptor transition with
  the globally optimal path. Structures with no shared anchor are invisible.
- The result is not the complete set of alternative isoforms and is not an
  exact second-best, third-best, or general k-best Viterbi result.
- Only one metadata-best reference/divergent hypothesis is retained per local
  cell. Alternatives needing a discarded metadata history can be missed.
- The bounded interval can truncate a distant transcript boundary.
- A locally optimal candidate may be rejected if it reaches a neighboring
  global locus; the search does not merge neighboring genes.
- Score deltas depend on the chosen equal genomic comparison interval.
- The current UniAnn HMM is forward-strand only. Alternative mode does not
  invent reverse-complement prediction, and no reverse-strand claim is made.
- A donor anchor fixes only that donor transition and an acceptor anchor only
  that acceptor transition. It does not fix both ends of the original intron.

These candidates are anchored exploratory hypotheses, not a complete isoform
catalog.

## Experimental gene-local k-best mode

The gene-local mode enumerates several complete transcript paths without
requiring a splice anchor from the global path:

```bash
scripts/uniann.sh -f sequence.fa -s scores.txt -p psauron.csv \
  --local-k-best \
  --local-k-output prediction.local_k_best.gff \
  --local-k-report prediction.local_k_best.tsv \
  --local-k-start-ratio 0.5
```

The default and current experiment use `K=5`, including the highest-scoring
path. For each global-best transcript, the fixed local interval begins at the
previous transcript's stop cell and ends immediately before the next
transcript's start cell; chromosome ends are used for the first and last
transcripts. This makes the neighboring global-best gene the hard boundary.

The reference ATG and upstream ATGs in the interval are eligible starts. An
upstream ATG is retained when its original start likelihood is at least the
configured fraction of the reference start likelihood; the ratio must satisfy
`0 < ratio <= 1`. The wrapper converts that ratio into a log-score drop using
the same dynamic `FACTOR` used to build
`out.atg.txt`. All eligible starts then receive the reference ATG transition
score inside the local DP. Only that score is inherited: genomic frame and
path metadata are recomputed at the candidate coordinate.

Paths have exactly three phases: pre-start `N`, coding `E*/I*`, and absorbing
post-stop `N`. A path cannot restart after a legal stop. Every retained path is
scored across the same complete local interval, including its `N` prefix and
suffix, so different start and stop coordinates remain comparable. Paths that
have not stopped before the next global-best gene are incomplete and excluded.

The probability recurrence follows Brown and Golod (2010), *Decoding HMMs
using the k best paths: algorithms and applications*: each destination cell
selects the K greatest extensions from the sorted predecessor hypotheses. The
paper's online compressed backpointer tree targets large K and long global
sequences. This implementation instead uses the paper's naive backpointer
matrix inside one local interval at a time, with active score/metadata
frontiers and a flat compact backpointer buffer.

UniAnn is not a plain first-order HMM: transition legality also depends on
`intron_len`, `exon_len`, `inter_len`, `exon_from`, and the predecessor state.
Each retained hypothesis therefore carries its own `PathMetadata`. The current
implementation still retains only K hypotheses per HMM state, not per fully
augmented metadata state. Consequently it is the direct practical K-best
