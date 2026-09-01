#!/usr/bin/env python3
"""Exact-CDS recovery of UniAnn local k-best paths against a reference annotation.

Reference : a GFF3/GTF with locus, transcript and CDS records, forward strand
            (for example the output of `gffread -M` on a curated annotation)
Query     : the GFF3 written by `uniann --local-k-best K --local-k-output`
Match     : identical ordered tuple of (start, end, phase) CDS segments

Every filter applied here is counted and reported; nothing is dropped silently.
"""

import argparse
import collections
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"           # inputs; gitignored by the repository
OUT = DATA / "eval"            # outputs; gitignored with the rest of data/


def parse_attributes(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_reference(path):
    """Returns loci, transcript->locus, transcript->CDS structure."""
    loci, tx_locus, cds = {}, {}, collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = parse_attributes(f[8])
        if f[2] == "locus":
            loci[a["ID"]] = (int(f[3]), int(f[4]))
        elif f[2] == "transcript":
            tx_locus[a["ID"]] = a.get("locus")
        elif f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]), int(f[4]), f[7]))
    structures = {t: tuple(sorted(v)) for t, v in cds.items()}
    return loci, tx_locus, structures


def read_kbest(path):
    """Returns locus->(start, end) and locus->[(rank, gene_id, structure)].

    Output layout: locus{I} > gene locus{I}.k{rank}.g{j} > transcript ...t1.
    A rank may contain several genes; each is returned as its own entry.
    """
    loci, tx_gene, cds = {}, {}, collections.defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = parse_attributes(f[8])
        if f[2] == "locus":
            loci[a["ID"]] = (int(f[3]), int(f[4]))
        elif f[2] == "transcript":
            tx_gene[a["ID"]] = a["Parent"]
        elif f[2] == "CDS":
            cds[a["Parent"]].append((int(f[3]), int(f[4]), f[7]))
    paths = collections.defaultdict(list)
    for t, segments in cds.items():
        gene = tx_gene[t]
        locus_id, rest = gene.split(".k", 1)
        rank = int(rest.split(".g")[0])
        paths[locus_id].append((rank, gene, tuple(sorted(segments))))
    for locus_id in paths:
        paths[locus_id].sort()
    return loci, paths


def intron_chain(structure):
    s = sorted(structure)
    return tuple((s[i][1], s[i + 1][0]) for i in range(len(s) - 1))


def map_loci_to_genes(loci, genes):
    """Each reference locus is mapped to the k-best search interval it overlaps most."""
    ordered = sorted(genes.items(), key=lambda kv: kv[1])
    mapping, ambiguous = {}, {}
    for locus_id, (ls, le) in loci.items():
        hits = [(min(le, ge) - max(ls, gs) + 1, g)
                for g, (gs, ge) in ordered if gs <= le and ls <= ge]
        if not hits:
            continue
        hits.sort(reverse=True)
        mapping[locus_id] = hits[0][1]
        if len(hits) > 1:
            ambiguous[locus_id] = [g for _, g in hits]
    return mapping, ambiguous


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", type=Path,
                   default=DATA / "NC_004354.4.p.CDS.uniq.gtf")
    p.add_argument("--kbest", type=Path,
                   default=DATA / "chrX_local_k_best.k20.gff")
    p.add_argument("--canonical-ids", type=Path, default=None,
                   help="optional list of reference transcript IDs surviving gffread -N")
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--max-k", type=int, default=0,
                   help="highest rank to consider; 0 uses every rank present")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    loci, tx_locus, ref = read_reference(args.reference)
    kbest_loci, paths = read_kbest(args.kbest)
    if not args.max_k:
        args.max_k = max((rank for entries in paths.values()
                          for rank, _g, _s in entries), default=1)
    mapping, ambiguous = map_loci_to_genes(loci, kbest_loci)

    canonical = None
    if args.canonical_ids and args.canonical_ids.exists():
        canonical = {l.strip() for l in open(args.canonical_ids) if l.strip()}

    structures_per_locus = collections.defaultdict(set)
    for t, s in ref.items():
        structures_per_locus[tx_locus[t]].add(s)

    # A reference structure may be predicted as the second gene of a path whose
    # interval belongs to a neighbouring locus, so candidates are indexed by
    # structure across all intervals rather than looked up per mapped locus.
    # CDS structures are genomic coordinate tuples, so they are unambiguous.
    rank_of_structure, rank_of_chain = {}, {}
    for entries in paths.values():
        for rank, _gene, structure in entries:
            if rank > args.max_k:
                continue
            if structure not in rank_of_structure or rank < rank_of_structure[structure]:
                rank_of_structure[structure] = rank
            chain = intron_chain(structure)
            if chain not in rank_of_chain or rank < rank_of_chain[chain]:
                rank_of_chain[chain] = rank

    rows = []
    for t, structure in sorted(ref.items()):
        locus_id = tx_locus[t]
        gene = mapping.get(locus_id)
        candidates = paths.get(gene, []) if gene else []
        exact_rank = rank_of_structure.get(structure)
        chain_rank = rank_of_chain.get(intron_chain(structure))
        rows.append({
            "transcript_id": t,
            "locus_id": locus_id,
            "kbest_gene": gene or "",
            "locus_mapped": int(gene is not None),
            "locus_ambiguous": int(locus_id in ambiguous),
            "exon_count": len(structure),
            "cds_length": sum(e - s + 1 for s, e, _ in structure),
            "structures_in_locus": len(structures_per_locus[locus_id]),
            "multi_structure_locus": int(len(structures_per_locus[locus_id]) > 1),
            "canonical_splice": "" if canonical is None else int(t in canonical),
            "candidates_in_locus": len(candidates),
            "exact_rank": exact_rank or "",
            "exact_match": int(exact_rank is not None),
            "intron_chain_rank": chain_rank or "",
            "intron_chain_match": int(chain_rank is not None),
        })

    with (args.out / "recovery_transcripts.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def summarise(subset, label):
        n = len(subset)
        exact = sum(r["exact_match"] for r in subset)
        chain = sum(r["intron_chain_match"] for r in subset)
        mapped = sum(r["locus_mapped"] for r in subset)
        return {
            "subset": label,
            "reference_structures": n,
            "in_mapped_locus": mapped,
            "exact_cds_match": exact,
            "exact_pct": round(100 * exact / n, 2) if n else 0.0,
            "intron_chain_match": chain,
            "intron_chain_pct": round(100 * chain / n, 2) if n else 0.0,
        }

    summary = [
        summarise(rows, "all reference structures"),
        summarise([r for r in rows if r["multi_structure_locus"]], "multi-structure loci only"),
        summarise([r for r in rows if not r["multi_structure_locus"]], "single-structure loci only"),
        summarise([r for r in rows if r["exon_count"] == 1], "single-CDS structures"),
        summarise([r for r in rows if r["exon_count"] > 1], "multi-CDS structures"),
    ]
    if canonical is not None:
        summary.append(summarise([r for r in rows if r["canonical_splice"] == 1],
                                 "canonical splice sites (gffread -N)"))

    by_rank = collections.Counter(r["exact_rank"] for r in rows if r["exact_match"])
    cumulative, running = {}, 0
    for k in range(1, args.max_k + 1):
        running += by_rank.get(k, 0)
        cumulative[k] = running

    counts = {
        "reference_loci": len(loci),
        "reference_transcripts": len(ref),
        "reference_unique_structures": len(set(ref.values())),
        "kbest_loci": len(kbest_loci),
        "kbest_paths": sum(len(v) for v in paths.values()),
        "loci_mapped_to_kbest_locus": len(mapping),
        "loci_unmapped": len(loci) - len(mapping),
        "loci_overlapping_multiple_genes": len(ambiguous),
        "exact_by_rank": {k: by_rank.get(k, 0) for k in range(1, args.max_k + 1)},
        "exact_cumulative": cumulative,
    }
    (args.out / "recovery_summary.json").write_text(
        json.dumps({"counts": counts, "summary": summary}, indent=2) + "\n")

    print(json.dumps(counts, indent=2))
    print()
    header = ["subset", "reference_structures", "in_mapped_locus",
              "exact_cds_match", "exact_pct", "intron_chain_match", "intron_chain_pct"]
    print("\t".join(header))
    for s in summary:
        print("\t".join(str(s[h]) for h in header))


if __name__ == "__main__":
    main()
