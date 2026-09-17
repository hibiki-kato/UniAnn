# Metadata-aware best Viterbi

`--metadata-state-viterbi` keeps the best score for each combination of coarse
state and transition-relevant metadata. It also replaces the full DP and
backpointer matrices with two score frontiers and a shared traceback tree.
The option is off by default. This change implements one best path; it does not
add k-best enumeration.

```sh
make -C src
src/uniann sequence.fa out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt \
  --metadata-state-viterbi > prediction.gff3
```

The wrapper `scripts/uniann.sh` accepts and forwards the same option. Input score
formats and precision are unchanged. The metadata decoder does not emit a DP
matrix log: its hidden state is a (coarse state, metadata) pair, so a per-base
matrix would have hundreds of columns and no reader. For the default decoder,
`--no-dp-dump` suppresses that log without changing decoding or GFF output; the
wrapper's `-n` now forwards it, so the matrix is not formatted only to be
discarded.

## Why the state includes metadata

The seven coarse states are N, E0–E2 and I0–I2. Transition rules also inspect
intron length, exon length and origin, and intergenic length and predecessor.
Two paths in I0 can therefore have different legal continuations: a recently
started intron can have the higher score while an older intron is the only one
long enough to use the next acceptor. Keeping only the higher-scoring I0 path
can discard the eventual optimum.

`tests/test_metadata_states.cpp` constructs this case and checks that both
histories survive until the acceptor. Its optimum, 950, agrees with a reference
that retains every distinct raw metadata key.

Only metadata distinctions that can affect a future transition are retained:

| Coarse state | Canonical metadata | Conservative reachable keys |
| --- | --- | ---: |
| N | Intergenic length capped at 28; whether predecessor is N | 28 |
| Each E state, entered from N | Exon length 3–103, capped at 103 | 101 |
| Each E state, entered from an intron | Exon length 1–5, capped at 5 | 5 |
| Each I state | Length 1–37, then three residue classes from 38 | 40 |

These thresholds follow the existing transition rules, including their
coordinate offsets and strict versus inclusive length comparisons. For example,
the single-exon stop rule checks `exon_len - 2 > 100`; lengths above 103 are
indistinguishable. Long introns still retain length modulo three.

There are at most `28 + 3 * (101 + 5) + 3 * 40 = 466` reachable keys under these
rules. The arithmetic ID mapping reserves **511 slots**, including unused
length/predecessor combinations. The decoder keeps one score and one traceback
head per reachable key, not one per possible complete path. Scores remain
`double`; frontier parent IDs use `uint16_t`. Ascending parent IDs and strict
score improvement give deterministic ties. Equal-score paths may differ from
the legacy decoder because its state space and tie order differ.

Few of those slots are live at once: on Dmel chrX the frontier holds 20 keys on
average and 73 at most. Each step therefore iterates an ascending list of the
active keys rather than all 511 slots, in both the score update and the
traceback tree, and resets only the entries it wrote. Sequence motifs at a
position (GT, AG, ATG, stop codon) are tested once per position rather than
once per candidate edge. These changes keep the tie order and produce the same
GFF bytes as the dense loops they replace.

## What remains in memory

The traceback tree stores parent links as `uint32_t` node IDs. It shares common
ancestry among surviving candidates and prunes a branch once neither a candidate
head nor a child needs it. New heads are acquired before old heads are released.
Consecutive equal states extend an unshared leaf instead of adding one node per
base. Common prefixes are periodically emitted into the final run list and
removed from the tree.

A single previous-state ID per current candidate would not recover a complete
path: the referenced earlier nodes must also survive. The tree retains exactly
that ancestry, with run compression and reuse of freed slots. It does not keep
all historical scores. Once the winner is known, replay calculates only the
scores at run boundaries needed by GFF output and verifies the total score.
The GFF writer preserves upstream coordinates, interval-average scores and the
special handling of the final segment.

Traceback memory depends on unresolved branching and run lengths. It is not
bounded by the frontier size alone, and long unresolved histories can still be
large. The node pool keeps its allocated capacity for reuse; final runs also
remain in memory until output. Node-ID exhaustion produces an error.

Inputs remain dense: five float emissions and one double site score per base
plus sequence storage, 29 bytes per base. On a long chromosome these are almost
all of the metadata decoder's memory (human chr1: 6.8 GiB peak, of which 7.2 GB
is input). Chunked input is a separate possible change. The current CLI expects
a single FASTA record for decoding, uses the forward strand of that input, and
requires sequence length to fit a signed `int`. Reverse-strand evaluation requires the corresponding transformed
inputs and reference, as in the existing workflow.

## Reading and testing the implementation

The files separate four responsibilities:

- `src/viterbi_model.h`: shared scoring rules and raw metadata updates.
- `src/metadata_states.h`: canonical keys, ID mapping and one score-update step.
- `src/shared_trace.h`: node ownership, branch pruning and run compression.
- `src/metadata_viterbi.h`: decoder orchestration and chosen-path score replay.

`src/uniann.cpp` selects the decoder and feeds state changes to the shared GFF
writer. The default decoder retains its historical DP/BT behavior, including
its finite unreachable-score sentinel. The new decoder uses negative infinity
for unreachable frontier entries and follows only allowed edges.

```sh
make -C src check
```

The checks cover 19 upstream GFF and full DP/BT golden outputs; randomized and
length-boundary comparisons against an uncanonicalized metadata reference;
full dense-parent traceback comparisons with and without run compression;
shared-tree comparisons against copied complete paths; and CLI and wrapper
behavior. Shared traceback and metadata decoding also passed Address
Sanitizer and UndefinedBehaviorSanitizer during PR preparation. These checks
validate the refactor against upstream and the new search/traceback against
separate reference implementations; they do not imply all inputs are covered.

## Measurements

Representative runs on 2026-09-17, compiled with the repository's `-O3 -march=native`
flags. RSS is process peak RSS from GNU `time -v`, converted from KiB to GiB.
Both Dmel modes used the same PR binary and prepared score inputs. The default
mode disabled the DP dump. Preprocessing and postprocessing are outside the
reported decoder measurements. Wall time includes reading the 1 GB emission
file, about 3.5 s of each run.

| Dmel chrX, 23,542,271 bases | SN (%) | PR (%) | F1 (%) | Peak RSS (GiB) | Wall time (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Default | 41.0544 | 69.1871 | 51.5312 | 4.59 | 6.07 |
| Metadata | 41.7274 | 69.8592 | 52.2472 | 0.64 | 9.08 |

Before the input, parsing and active-key changes (runs on 2026-09-08 with the
same accuracy counts), the same inputs took 19.70 s / 8.12 GiB in default mode
and 41.36 s / 1.34 GiB in metadata mode. Both decoders' GFF outputs are
byte-identical before and after those changes.

The reference has 1,783 distinct forward-strand CDS structures. Default decoding
produced 1,058 predictions with 732 exact matches; metadata decoding produced
1,065 with 744 matches. Matching requires sequence, strand, all CDS coordinates
and phases to agree. The existing wrapper's short one-/two-exon filter was
applied equally. `gffcompare --strict-match -e 0` on the same outputs gives
transcript-level SN/PR 41.1/69.2 (default) and 41.7/69.9 (metadata); on the
first 100 Mb of human chr1, 10.9/52.2 and 11.2/53.2. This is one scoring
configuration, not a claim of a universal accuracy or speed gain.

Human chr1 (`NC_000001.11`, 248,956,422 bases) completed with the production
metadata decoder at **6.78 GiB** peak RSS (7,111,040 KiB), **97.4 seconds**
(438.67 seconds and 14.20 GiB before the changes above). Its GFF,
including scores, was byte-identical between the two versions and to the result
of rendering an independently checkpoint-restored metadata path with the
unmodified upstream GFF writer. Dmel passed the same full-output comparison.
The human best-path score was `85168357.143762767`. The default decoder was
stopped by a 48 GiB address-space limit on the same input, so a paired human
accuracy evaluation is not reported here.

The earlier checkpoint experiment and its large local inputs are not bundled
with this PR. The self-contained small reference tests above are included.
