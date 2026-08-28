#!/usr/bin/env python3
"""Small deterministic integration tests for UniAnn's legacy and alt modes."""

import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "src" / "uniann"
LEGACY_STDOUT_SHA256 = "55a78a0506c5472da059fbf0b28b529ee3646747cf5f8c1fd11d1585798d7a5d"
LEGACY_STDERR_SHA256 = "229d9600474847f410a43e6ff583396ea49415b46594f36f46054996e6fe7f46"


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def write_inputs(directory, name, sequence, emissions, gt=(), ag=(), atg=(), stop=()):
    directory = Path(directory)
    fasta = directory / f"{name}.fa"
    fasta.write_text(f">{name}\n{sequence}\n")
    emit = directory / f"{name}.emit"
    with emit.open("w") as handle:
        for position, values in enumerate(emissions):
            handle.write(str(position) + "\t" + "\t".join(map(str, values)) + "\n")
    paths = []
    for suffix, scores in (("gt", gt), ("ag", ag), ("atg", atg), ("stop", stop)):
        path = directory / f"{name}.{suffix}"
        path.write_text("".join(f"{position}\t{score}\n" for position, score in scores))
        paths.append(path)
    return [fasta, emit, *paths]


def run(inputs, *options):
    return subprocess.run(
        [str(BINARY), *(str(path) for path in inputs), *map(str, options)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def make_spliced_fixture(directory, kind="donor"):
    length = 240
    sequence = ["A"] * length
    motifs = [(0, "ATG"), (60, "GT"), (90, "AG"), (200, "TAA")]
    gt = [(60, 100)]
    ag = [(90, 100)]
    if kind == "donor":
        motifs.append((51, "GT"))
        gt.insert(0, (51, 95))
    elif kind == "acceptor":
        motifs.append((99, "AG"))
        ag.append((99, 95))
    for position, motif in motifs:
        sequence[position : position + len(motif)] = motif

    emissions = []
    for position in range(length):
        values = [0.0] + [-1000.0] * 6
        if 2 <= position <= 60:
            values[1] = 1.0  # E0
        if 61 <= position <= 99:
            values[4] = 1.0  # I0
        if kind == "donor" and 52 <= position <= 60:
            values[4] = 1.0
        if 91 <= position <= 201:
            values[3] = 1.0  # E2 after a 32- or 41-base modeled intron
        emissions.append(values)
    return write_inputs(
        directory,
        f"alternative_{kind}",
        "".join(sequence),
        emissions,
        gt=gt,
        ag=ag,
        atg=[(0, 100)],
        stop=[(200, 100)],
    )


def make_three_exon_fixture(directory, mode):
    length = 280
    sequence = ["A"] * length
    motifs = [(0, "ATG"), (60, "GT"), (90, "AG"),
              (150, "GT"), (180, "AG"), (220, "TAA")]
    if mode == "duplicate":
        motifs.append((51, "GT"))
    else:
        motifs.append((230, "TAA"))
    for position, motif in motifs:
        sequence[position : position + len(motif)] = motif

    emissions = []
    for position in range(length):
        values = [0.0] + [-1000.0] * 6
        if 2 <= position <= 60:
            values[1] = 1.0
        if (52 if mode == "duplicate" else 61) <= position <= 180:
            values[4] = 1.0 if mode in {"duplicate", "inclusion"} else 0.0
        if 91 <= position <= 150:
            values[3] = 1.0 if mode != "inclusion" else 0.0
        if 151 <= position <= 180:
            values[6] = 1.0 if mode != "inclusion" else 0.0
        if 181 <= position <= 221:
            values[2] = 1.0 if mode != "inclusion" else 0.0
        if mode in {"skipping", "inclusion"} and 181 <= position <= 231:
            values[3] = 1.0
        emissions.append(values)

    if mode == "duplicate":
        gt = [(51, 95), (60, 100), (150, 100)]
        ag = [(90, 100), (180, 100)]
        stop = [(220, 100)]
    elif mode == "skipping":
        gt = [(60, 100), (150, 100)]
        ag = [(90, 100), (180, 100)]
        stop = [(220, 100), (230, 90)]
    elif mode == "inclusion":
        gt = [(60, 100), (150, -50)]
        ag = [(90, -50), (180, 100)]
        stop = [(220, 90), (230, 100)]
    else:
        raise ValueError(mode)
    return write_inputs(directory, mode, "".join(sequence), emissions,
                        gt=gt, ag=ag, atg=[(0, 100)], stop=stop)


def parse_gff(path):
    features = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == 9, line
        attributes = {}
        for item in fields[8].split(";"):
            key, value = item.split("=", 1)
            attributes[key] = value
        features.append((fields, attributes))
    return features


def assert_valid_gff(path):
    features = parse_gff(path)
    ids = set()
    parents = []
    for fields, attributes in features:
        start, end = int(fields[3]), int(fields[4])
        assert 1 <= start <= end
        assert fields[6] == "+"
        if fields[2] == "CDS":
            assert fields[7] in {"0", "1", "2"}
        else:
            assert fields[7] == "."
        if "ID" in attributes:
            assert attributes["ID"] not in ids
            ids.add(attributes["ID"])
        if "Parent" in attributes:
            parents.append(attributes["Parent"])
    assert parents and all(parent in ids for parent in parents)
    return features


def test_legacy_golden():
    with tempfile.TemporaryDirectory() as temporary:
        sequence = "ATGCGTACGTACGTACGTAA"
        emissions = [[0] * 7 for _ in sequence]
        inputs = write_inputs(temporary, "tiny", sequence, emissions)
        completed = run(inputs)
        assert completed.returncode == 0
        assert sha256(completed.stdout) == LEGACY_STDOUT_SHA256
        assert sha256(completed.stderr) == LEGACY_STDERR_SHA256


def run_alt_fixture(kind):
    temporary = tempfile.TemporaryDirectory()
    directory = Path(temporary.name)
    inputs = make_spliced_fixture(directory, kind)
    gff = directory / "alternative.gff3"
    report = directory / "alternative.tsv"
    completed = run(
        inputs,
        "--alternative-splicing",
        "--alt-output",
        gff,
        "--alt-report",
        report,
    )
    assert completed.returncode == 0, completed.stderr.decode()
    rows = list(csv.DictReader(report.open(), delimiter="\t"))
    return temporary, completed, gff, rows


def test_alternative_donor_and_format():
    temporary, completed, gff, rows = run_alt_fixture("donor")
    try:
        kept = [row for row in rows if row["kept_after_deduplication"] == "1"]
        assert len(kept) == 1
        assert kept[0]["classification"] == "alternative_donor"
        assert kept[0]["intron_chain"] == "51-90"
        assert float(kept[0]["score_delta"]) == -5.0
        assert (float(kept[0]["candidate_total_score"]) -
                float(kept[0]["reference_interval_score"]) ==
                float(kept[0]["score_delta"]))
        assert kept[0]["status"] == "complete" and not kept[0]["failure_reason"]
        features = assert_valid_gff(gff)
        transcripts = [a["ID"] for f, a in features if f[2] == "transcript"]
        assert transcripts == ["UniAnnGene000001.t1", "UniAnnGene000001.t2"]
        assert b"Number of alternative donors: 1\n" in completed.stderr
    finally:
        temporary.cleanup()


def test_alternative_acceptor():
    temporary, completed, gff, rows = run_alt_fixture("acceptor")
    try:
        kept = [row for row in rows if row["kept_after_deduplication"] == "1"]
        assert len(kept) == 1
        assert kept[0]["classification"] == "alternative_acceptor"
        assert kept[0]["intron_chain"] == "60-99"
        assert_valid_gff(gff)
        assert b"Number of alternative acceptors: 1\n" in completed.stderr
    finally:
        temporary.cleanup()


def test_one_unambiguous_transcript():
    with tempfile.TemporaryDirectory() as temporary:
        length = 170
        sequence = ["A"] * length
        sequence[0:3] = "ATG"
        sequence[120:123] = "TAA"
        emissions = []
        for position in range(length):
            values = [0.0] + [-1000.0] * 6
            if 2 <= position <= 121:
                values[1] = 1.0
            emissions.append(values)
        inputs = write_inputs(temporary, "single", "".join(sequence), emissions,
                              atg=[(0, 100)], stop=[(120, 100)])
        gff = Path(temporary) / "single.gff3"
        report = Path(temporary) / "single.tsv"
        completed = run(inputs, "--alternative-splicing", "--alt-output", gff,
                        "--alt-report", report)
        assert completed.returncode == 0
        assert b"Number of splice-site anchors: 0\n" in completed.stderr
        features = assert_valid_gff(gff)
        assert sum(f[2] == "transcript" for f, _ in features) == 1
        assert len(list(csv.DictReader(report.open(), delimiter="\t"))) == 0


def test_local_cutoff_is_incomplete_and_not_primary():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        inputs = make_spliced_fixture(directory, "donor")
        gff = directory / "cutoff.gff3"
        report = directory / "cutoff.tsv"
        completed = run(inputs, "--alternative-splicing", "--alt-output", gff,
                        "--alt-report", report, "--alt-include-incomplete",
                        "--alt-max-distance", "10")
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="\t"))
        assert rows and all(row["status"] == "maximum_distance_reached" for row in rows)
        features = assert_valid_gff(gff)
        assert sum(f[2] == "transcript" for f, _ in features) == 1


def test_duplicate_discovery_is_deduplicated():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        inputs = make_three_exon_fixture(directory, "duplicate")
        gff = directory / "duplicate.gff3"
        report = directory / "duplicate.tsv"
        completed = run(inputs, "--alternative-splicing", "--alt-output", gff,
                        "--alt-report", report)
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="\t"))
        alternatives = [r for r in rows if r["classification"] == "alternative_donor"]
        assert len(alternatives) >= 2
        assert sum(r["kept_after_deduplication"] == "1" for r in alternatives) == 1
        features = assert_valid_gff(gff)
        assert sum(f[2] == "transcript" for f, _ in features) == 2


def test_exon_skipping_and_intron_retention():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        inputs = make_three_exon_fixture(directory, "skipping")
        gff = directory / "skipping.gff3"
        report = directory / "skipping.tsv"
        completed = run(inputs, "--alternative-splicing", "--alt-output", gff,
                        "--alt-report", report)
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="\t"))
        classes = {r["classification"] for r in rows
                   if r["kept_after_deduplication"] == "1"}
        assert "exon_skipping" in classes
        assert "intron_retention" in classes
        assert all(r["reference_interval_score"] != "." for r in rows)
        assert_valid_gff(gff)


def test_alternative_exon_inclusion():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        inputs = make_three_exon_fixture(directory, "inclusion")
        gff = directory / "inclusion.gff3"
        report = directory / "inclusion.tsv"
        completed = run(inputs, "--alternative-splicing", "--alt-output", gff,
                        "--alt-report", report)
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="\t"))
        assert any(r["classification"] == "alternative_exon" and
                   r["kept_after_deduplication"] == "1" for r in rows)
        assert_valid_gff(gff)


def test_cli_rejects_alt_options_without_mode():
    with tempfile.TemporaryDirectory() as temporary:
        sequence = "A" * 20
        inputs = write_inputs(temporary, "cli", sequence, [[0] * 7 for _ in sequence])
        completed = run(inputs, "--alt-output", Path(temporary) / "bad.gff")
        assert completed.returncode == 1
        assert b"require --alternative-splicing" in completed.stderr



def run_local_kbest_fixture(kind):
    temporary = tempfile.TemporaryDirectory()
    directory = Path(temporary.name)
    inputs = make_spliced_fixture(directory, kind)
    gff = directory / "local_kbest.gff3"
    report = directory / "local_kbest.tsv"
    completed = run(
        inputs,
        "--local-k-best",
        "--local-k-output",
        gff,
        "--local-k-report",
        report,
    )
    assert completed.returncode == 0, completed.stderr.decode()
    rows = list(csv.DictReader(report.open(), delimiter="	"))
    return temporary, completed, gff, rows

def test_local_kbest_alternate_donor_and_ranking():
    temporary, completed, gff, rows = run_local_kbest_fixture("donor")
    try:
        assert len(rows) > 0
        classes = [r["classification"] for r in rows]
        assert "alternative_donor" in classes
        ranks = [int(r["rank"]) for r in rows]
        assert ranks == list(range(1, len(rows) + 1))
        features = assert_valid_gff(gff)
        transcripts = [a["ID"] for f, a in features if f[2] == "transcript"]
        assert len(transcripts) >= 2
    finally:
        temporary.cleanup()

def test_local_kbest_cli_validation():
    with tempfile.TemporaryDirectory() as temporary:
        sequence = "A" * 20
        inputs = write_inputs(temporary, "cli", sequence, [[0] * 7 for _ in sequence])
        completed = run(inputs, "--local-k-output", Path(temporary) / "bad.gff")
        assert completed.returncode == 1

        completed = run(inputs, "--local-k-best",
                        "--local-k-start-score-drop", "-1",
                        "--no-dp-dump")
        assert completed.returncode == 1
def test_local_kbest_unique_dedup():
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        inputs = make_three_exon_fixture(directory, "duplicate")
        gff = directory / "local_kbest.gff3"
        report = directory / "local_kbest.tsv"
        completed = run(inputs, "--local-k-best", "10", "--local-k-output", gff, "--local-k-report", report)
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="	"))
        chains = [r["intron_chain"] for r in rows]
        assert len(chains) == len(set(chains)), "Paths must be unique"
        features = assert_valid_gff(gff)
        transcripts = [a["ID"] for f, a in features if f[2] == "transcript"]
        assert len(transcripts) <= 10

def test_local_kbest_candidate_start_inherited_score():
    with tempfile.TemporaryDirectory() as temporary:
        length = 170
        sequence = ["A"] * length
        sequence[0:3] = "ATG"
        sequence[10:13] = "ATG"
        sequence[120:123] = "TAA"
        emissions = []
        for position in range(length):
            values = [0.0] + [-1000.0] * 6
            if 2 <= position <= 121:
                values[1] = 1.0
            emissions.append(values)
        inputs = write_inputs(temporary, "single", "".join(sequence), emissions,
                              atg=[(0, 100), (10, 80)], stop=[(120, 100)])
        gff = Path(temporary) / "single.gff3"
        report = Path(temporary) / "single.tsv"
        completed = run(inputs, "--local-k-best", "--local-k-output", gff,
                        "--local-k-report", report, "--local-k-start-score-drop", "50")
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="	"))
        inherited = [float(r["inherited_start_score"]) for r in rows if r["actual_start_score"] != "."]
        actual = [float(r["actual_start_score"]) for r in rows if r["actual_start_score"] != "."]
        assert any(abs(i - a) > 10 for i, a in zip(inherited, actual)), "Inherited score should differ from actual for alternates"

def test_local_kbest_stopping_before_next_gene():
    with tempfile.TemporaryDirectory() as temporary:
        length = 300
        sequence = ["A"] * length
        sequence[0:3] = "ATG"
        sequence[120:123] = "TAA"
        sequence[150:153] = "ATG"
        sequence[270:273] = "TAA"
        emissions = []
        for position in range(length):
            values = [0.0] + [-1000.0] * 6
            if 2 <= position <= 121 or 152 <= position <= 271:
                values[1] = 1.0
            emissions.append(values)
        inputs = write_inputs(temporary, "multi", "".join(sequence), emissions,
                              atg=[(0, 100), (150, 100)], stop=[(120, 100), (270, 100)])
        gff = Path(temporary) / "multi.gff3"
        report = Path(temporary) / "multi.tsv"
        completed = run(inputs, "--local-k-best", "--local-k-output", gff,
                        "--local-k-report", report)
        assert completed.returncode == 0
        rows = list(csv.DictReader(report.open(), delimiter="\t"))
        first_gene = [row for row in rows
                      if row["reference_gene_index"] == "1"]
        assert first_gene
        reference_start = int(first_gene[0]["reference_start_1based"])
        assert all(int(row["candidate_end_1based"]) >= reference_start
                   for row in first_gene)
        features = assert_valid_gff(gff)
        transcripts = [f for f, _ in features if f[2] == "transcript"]
        ends = [f[4] for f in transcripts]
        assert all(int(end) <= length for end in ends)

def main():
    tests = [value for name, value in sorted(globals().items()) if name.startswith("test_")]
    for test in tests:
        test()
        print(f"PASS {test.__name__}")
    print(f"PASS {len(tests)} integration tests")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise
