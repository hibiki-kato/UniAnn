# UniAnn alternative-splicing design note

## Existing implementation

`src/uniann.cpp` contains the complete predictor used by `scripts/uniann.sh`.
The wrapper generates six score files, invokes the C++ executable, redirects its
stderr DP dump to `out.err`, and pipes its stdout through `gffread -F`.

The HMM has exactly seven forward-strand states: `N`, `E0`, `E1`, `E2`, `I0`,
`I1`, and `I2`. `N` is intergenic/noncoding. `E*` is coding exon sequence and
`I*` is intron sequence; the suffix records coding frame. There are no reverse-
strand, UTR, explicit start/end, or special first/final exon states.

The DP is initialized by `init_dp()` at zero-based sequence position 0, with
only `N` reachable. `run_viterbi()` advances left-to-right and stores one
`DPCell` per position/state. `viterbi_termination()` selects the greatest score
among all states at the last sequence position, retaining the original strict
`>` state-order tie break. `viterbi_backtrace()` follows the stored predecessor
state from that terminal cell.

Every cell stores cumulative score, predecessor state, current intron length,
current exon length, current intergenic length, and the state from which the
current exon began. Lengths and `exon_from` affect later legality, so a state
label alone is not a sufficient local-search condition.

All sequence, score-array, and DP indices are zero-based. A sparse GT, AG, ATG,
or stop score is indexed by the first motif base. An exon-to-intron transition
entering DP position `i` consumes the GT beginning at `i-1`; an intron-to-exon
transition consumes the AG beginning at `i-1`; start and stop transitions
entering position `i` consume the codon beginning at `i-2`. Named coordinate
helpers will preserve these conventions.

Self transitions have score zero. Destination emissions are added at every DP
position. Cross-state transitions are normally impossible and become legal only
through the GT, AG, ATG, or stop rules. The minimum intron, exon, intergenic,
and single-exon lengths are 30, 30, 30, and 100 bases. First exons are exempted
from the internal minimum-exon check. Intron exit uses intron length modulo
three plus the source/destination frame states. Start and stop codon coordinates
must agree with the exon frame. An in-frame stop blocks the exon self transition.

The raw writer calls its stream GFF3 but emits only `CDS` and `intron` features,
always on `+`, with phase `.`, zero-based-derived coordinates, and a `Parent`
attribute whose parent feature is not emitted. It has no exon/transcript/gene objects
and cannot represent multiple transcripts under one gene. The wrapper's
`gffread` pass is therefore part of the installed output behavior. Alternative
mode will introduce a reusable path-to-transcript representation and a valid,
deterministic GFF3 hierarchy while leaving the legacy stdout path byte-for-byte
unchanged when the mode is absent.

## Refactor and search design

The global loop, initialization, terminal selection, backtrace, strict tie
breaking, legacy stdout, and legacy stderr remain compatibility surfaces. The
biological recurrence will move behind one transition evaluator that returns
legality, additive score, updated metadata, and a failure reason. Global DP,
local DP, explicit rescoring, and validation will all call it.

After global backtrace, consecutive `E* -> I*` and `I* -> E*` transitions are
recorded as donor and acceptor anchors. The stored genomic coordinate is always
the zero-based first base of the GT or AG dinucleotide; the one-based report
coordinate is that value plus one. Adjacent donor/acceptor sites on the global
path are additionally paired as intron metadata. Each anchor retains its exact
DP transition position, states, frame, score, strand, and reference transcript.

The upstream implementation uses the allowed constrained-forward fallback.
Within `[anchor - max_distance, anchor]`, a virtual intergenic source represents
all legal transcript starts in one pass. The DP must end with the exact anchored
left-to-right transition. Downstream continuation begins with the resulting
right state and complete metadata, considers every legal transcript exit to
`N` through the interval end, and chooses the greatest complete endpoint rather
than the first endpoint.

Conditioning only on a transition already in a global optimum ordinarily
reconstructs that optimum by the Viterbi optimal-substructure principle. To make
the experimental search capable of returning an alternative without pretending
to enumerate k-best paths, local cells carry two deterministic hypotheses:
`reference_splice_structure` and `diverged_splice_structure`. Divergence is set
only by a splice transition difference relative to the associated global
transcript. The best complete divergent hypothesis is the candidate; the
reference hypothesis is retained for consistency checks. This is an explicit
constraint, not a score perturbation, and adds only a constant factor to local
DP complexity.

Upstream and downstream paths are joined in increasing genomic order with the
anchor transition represented once. A validator checks contiguity, bounds,
every transition and emission, motif coordinates, frame, codon and length rules,
boundary completeness, exact anchor presence, annotation conversion, and an
explicit score recomputation tolerance of `1e-8` (scale-aware when needed).

Paths convert to a shared `TranscriptAnnotation` containing exon, CDS, and
intron intervals. Full-structure and canonical ordered-intron-chain keys support
deterministic deduplication. Candidates are assigned to the anchor's reference
locus and rejected if they cross a neighboring reference gene. A candidate is
an alternative only when its intron chain differs; conservative classifications
are alternative donor, alternative acceptor, exon skipping, alternative exon,
intron retention, or complex.

Alternative mode writes an additional valid GFF3 file. Reference `.t1` and
score-ordered unique alternatives share a stable gene ID. The legacy normal
stdout remains unchanged. Optional TSV diagnostics may also retain invalid,
incomplete, identical, and deduplicated discoveries.

## Scope limitation

The search can find only structures sharing at least one selected donor or
acceptor transition with the global path. It cannot find an unrelated gene or
an isoform with no anchored splice site in common, is not exact k-best Viterbi,
and can miss structures beyond the bounded interval or structures requiring a
different metadata hypothesis discarded by the one-best-per-cell recurrence.
The current model and implementation are forward-strand only; no reverse-strand
correctness claim is made.

## Phase 1 checklist

1. Initialization: `init_dp()`, position zero, only `N` reachable.
2. Recurrence: `run_viterbi()`, positions `1..L-1`, all 49 state pairs.
3. Terminal state: maximum cell at `L-1`, strict `>` tie break in state order.
4. Backtrace: one predecessor state per cell, followed right-to-left.
5. States: exactly `N`, `E0..E2`, and `I0..I2` as described above.
6. Coordinates: FASTA, score arrays, DP, and internal intervals are zero-based.
7. Donor/acceptor: sparse score index is the first GT/AG base; the transition
   enters the DP cell on the motif's second base.
8. GT/AG: their sparse scores replace otherwise-forbidden E/I transitions.
9. Start/stop: sparse score index is the first codon base; the transition enters
   on the third base. Starts before position 25 retain the legacy score `1`.
10. Segment transitions: self score is zero, cross-state score comes only from
    the motif/boundary rules, and every destination emission is added.
11. Lengths: intron/exon/intergenic/single limits are 30/30/30/100; the first
    exon exemption and strict single-exon comparison are preserved.
12. Frame: encoded by E/I suffix, splice length modulo three, and genomic codon
    start modulo three.
13. Cell metadata: score, backpointer, intron/exon/intergenic lengths, and the
    state from which the current exon began.
14. Strand: literal `+`; no reverse-strand state graph exists.
15. Path output: consecutive labels become state runs with motif-specific end
    adjustments in `write_gff_from_path()`.
16. Format: raw custom GFF3-like output, normalized by wrapper `gffread -F`.
17. IDs: legacy `Parent=<FASTA-basename>.<global region counter>` only.
18. Features: legacy raw output emits CDS and intron, not exon/transcript/gene.
19. Phase: legacy raw output writes `.`, while `gffread` infers phase; the new
    hierarchical writer computes phase directly from cumulative CDS length.
20. Multiple transcripts: unsupported by the legacy writer; alternative mode
    uses the new transcript representation and hierarchical writer.

For reverse optimization, state/frame are reversible, segment lengths are
reconstructible only when the boundary is known, and `exon_from`, intergenic
history, and start/stop history cannot be recovered from a current state alone.
That classification is why the implementation uses constrained forward search.
