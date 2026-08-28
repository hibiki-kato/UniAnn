# Synthetic test fixtures

`tests/test_uniann.py` generates every fixture in a temporary directory. It does
not commit biological data or leave scratch score arrays behind.

- `tiny` is the 20-base sequence `ATGCGTACGTACGTACGTAA`; all seven emissions
  are zero and sparse arrays are empty. Its pristine stdout/stderr SHA-256
  hashes freeze the legacy annotation, full DP dump, and exit status.
- `single` is 170 bases with ATG at zero and TAA at 120. E0 emissions are high
  from DP 2 through 121. It has one complete, unambiguous, intronless transcript
  and therefore no anchors or alternative intron chains.
- `alternative_donor` is 240 bases with ATG 0, donors 51/60, acceptor 90, and
  stop 200. The donor scores are 95/100. Compatible E0, I0, and E2 emissions
  make donor 60 global and donor 51 the best acceptor-anchored alternative.
  Expected chains are `(60,90)` and `(51,90)`.
- `alternative_acceptor` uses donor 60, acceptors 90/99, and scores 100/95.
  The nine-base shift preserves frame. Expected chains are `(60,90)` and
  `(60,99)`.
- `duplicate` has global introns `(60,90),(150,180)` plus donor 51. The same
  `(51,90),(150,180)` alternative is recovered from multiple downstream anchors
  but only one representative is written.
- `skipping` has global exons separated by `(60,90),(150,180)`. Staying in I0
  from donor 60 to acceptor 180 yields the expected exon-skipping chain
  `(60,180)`. A second anchored path remains exonic across the second reference
  intron and exercises intron-retention classification.
- `inclusion` makes `(60,180)` globally optimal while retaining lower-scoring
  donor/acceptor evidence at 90/150. The divergent hypothesis recovers the
  three-exon `(60,90),(150,180)` structure.
- `cutoff` reuses the donor fixture with a 10-base search bound. Both anchors
  begin inside modeled transcript segments, so they are explicitly incomplete
  and absent from primary GFF3.

All unspecified sequence bases are `A`. Unspecified motif scores are negative
infinity. The intended state paths use E0 before the first donor, I0 through the
first intron, E2 after a 32- or 41-base modeled first intron, I2 through the
second intron, and E1 after its 32-base modeled length. Non-path emissions are
`-1000`; path emissions are normally `1` (or `0` where alternatives must be
compared only by splice score).

`tests/test_classification.cpp` adds direct, inspectable unit cases for:

- alternative donor and acceptor;
- exon skipping, alternative exon inclusion, intron retention, and complex;
- donor/acceptor/exon coordinate conversion at boundaries;
- multi-segment CDS phase;
- reading-frame rejection;
- missing or wrong-frame start/stop rejection;
- minimum exon, intron, intergenic, and single-exon rejection;
- authoritative transition-score use.

The Python GFF3 check verifies column count, coordinates, valid phases, unique
IDs, and that every child references an existing parent.
