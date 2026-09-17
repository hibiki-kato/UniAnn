#!/usr/bin/env python3
"""Paired global decoder benchmark using precomputed input scores.

Evaluates the exact CDS structure prediction performance (SN/PR/F1)
for a baseline and a metadata Viterbi model. Post-processes the decoded
GFF3 using gffread with a >200 bp filter for 1/2 exon predictions.
Maintains consistent hashes across runs to ensure comparability.

Example:
    python benchmarks/compare_viterbi.py \
        --data data/hsap \
        --out results/hsap_eval \
        --variants baseline metadata
"""

import argparse
import collections
import csv
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    """Computes the SHA-256 hash of a file."""
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def execute(command, stdout, stderr):
    """Executes a command, writing standard output and error to files."""
    with stdout.open('wb') as out, stderr.open('wb') as err:
        subprocess.run([str(x) for x in command], stdout=out, stderr=err,
                       env={**os.environ, 'LC_ALL': 'C'}, check=True)


def postprocess(raw, directory, gffread, min_cds):
    """Filters short 1/2 exon predictions identically to the production wrapper."""
    tlf = directory / 'all.tlf'
    execute([gffread, raw, '--tlf'], tlf, directory / 'gffread_tlf.stderr')

    filtered = directory / 'filtered.tlf'
    pattern = re.compile(r'exonCount=(1|2);exons=(\S+);CDS=(\d+):(\d+);CDSphase=\d')

    with filtered.open('w') as out:
        for line in tlf.read_text().splitlines(keepends=True):
            fields = line.split('\t')
            match = None
            if len(fields) > 8:
                match = pattern.search(fields[8])

            if match:
                length = 0
                for start, end in re.findall(r'(\d+)-(\d+)', match.group(2)):
                    length += int(end) - int(start) + 1

                if length <= min_cds:
                    continue

            out.write(line)

    annotation = directory / 'annotation.gff3'
    execute([gffread, filtered], annotation, directory / 'gffread.stderr')
    return annotation


def structures(path):
    """Reads CDS segments, grouped and deduped by transcript/parent."""
    cds = collections.defaultdict(list)

    for line in path.read_text().splitlines():
        if line.startswith('#'):
            continue

        fields = line.split('\t')
        if len(fields) < 9 or fields[2] != 'CDS':
            continue

        attributes_raw = fields[8]
        parents = []

        # Parse GFF3 Parent= or GTF transcript_id
        if 'Parent=' in attributes_raw:
            attributes = {}
            for item in attributes_raw.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attributes[key] = value
            if 'Parent' in attributes:
                parents = attributes['Parent'].split(',')

        elif 'transcript_id' in attributes_raw:
            match = re.search(r'transcript_id\s+"([^"]+)"', attributes_raw)
            if match:
                parents = [match.group(1)]

        if not parents or any(not parent for parent in parents):
            raise ValueError(f"Missing Parent or transcript_id in CDS row: {line}")

        seq = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        phase = fields[7]

        for parent in parents:
            cds[(seq, strand, parent)].append((start, end, phase))

    result = set()
    for (seq, strand, _), segments in cds.items():
        sorted_segments = tuple(sorted(segments))
        result.add((seq, strand, sorted_segments))

    return result


def exact_metrics(reference, query):
    """Computes precision, recall, and F1 on exact structure matches."""
    truth = structures(reference)
    pred = structures(query)
    tp = len(truth & pred)

    sn = 100 * tp / len(truth) if truth else 0.0
    pr = 100 * tp / len(pred) if pred else 0.0
    f1 = 200 * tp / (len(truth) + len(pred)) if (truth or pred) else 0.0

    return {
        'reference': len(truth),
        'predictions': len(pred),
        'tp': tp,
        'fp': len(pred - truth),
        'fn': len(truth - pred),
        'sn': sn,
        'pr': pr,
        'f1': f1
    }


def parse_time_v(log_path):
    """Parses wall clock seconds and peak RSS from `/usr/bin/time -v` output."""
    text = log_path.read_text()
    wall = None
    rss = None

    for line in text.splitlines():
        if 'Elapsed (wall clock) time' in line:
            val = line.split('): ')[1].strip()
            parts = val.split(':')
            if len(parts) == 2:
                wall = int(parts[0]) * 60 + float(parts[1])
            elif len(parts) == 3:
                wall = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif 'Maximum resident set size' in line:
            rss = int(line.split('): ')[1].strip())

    if wall is None or rss is None:
        raise ValueError(f"Missing resource measurements in {log_path}")
    return {"wall_seconds": wall, "max_rss_kib": rss}


def evaluate(args, variant, inputs):
    """Runs a single decoder variant and computes its evaluation metrics."""
    directory = args.out / variant
    directory.mkdir(exist_ok=True)
    # A failed rerun must not leave an older success in the comparison.
    (directory / 'metrics.json').unlink(missing_ok=True)

    binary = args.baseline_binary if variant == 'baseline' else args.metadata_binary
    command = [str(binary), *map(str, inputs), '--no-dp-dump']
    if variant == 'metadata':
        command.append('--metadata-state-viterbi')

    input_hashes = {}
    for p in [*inputs, args.reference]:
        input_hashes[str(p)] = digest(p)

    manifest = {
        'command': command,
        'binary_sha256': digest(binary),
        'inputs': input_hashes,
        'gffread_sha256': digest(args.gffread),
        'min_cds': args.min_cds
    }

    manifest_path = directory / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')

    print(f'{variant}: decoding', flush=True)
    time_log = directory / 'time.log'
    execute(['/usr/bin/time', '-v', '-o', time_log, *command],
            directory / 'raw.gff3', directory / 'decoder.stderr')

    annotation = postprocess(directory / 'raw.gff3', directory, args.gffread, args.min_cds)

    result = {
        'variant': variant,
        'exact_cds': exact_metrics(args.reference, annotation),
        'resources': parse_time_v(time_log)
    }

    metrics_path = directory / 'metrics.json'
    metrics_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


def summarize(args):
    """Collects individual run metrics and outputs a joint comparison summary."""
    results = []
    for variant in ('baseline', 'metadata'):
        metrics_file = args.out / variant / 'metrics.json'
        if metrics_file.exists():
            results.append(json.loads(metrics_file.read_text()))

    rows = []
    for result in results:
        variant = result['variant']
        m = result['exact_cds']
        rows.append({
            'variant': variant,
            'method': 'phase_aware_exact',
            'level': 'CDS_structure',
            'sn': m['sn'],
            'pr': m['pr'],
            'f1': m['f1']
        })

    summary = {
        'results': results,
        'note': 'SN, PR and F1 are percentages. Phase-aware exact CDS F1 uses integer structure counts. '
                'Reference must cover the same sequence and strand as the decoder input.'
    }

    if len(results) == 2:
        manifests = []
        for variant in ('baseline', 'metadata'):
            manifest_file = args.out / variant / 'manifest.json'
            manifests.append(json.loads(manifest_file.read_text()))

        for key in ('inputs', 'gffread_sha256', 'min_cds'):
            if manifests[0][key] != manifests[1][key]:
                raise ValueError(f'Paired evaluation differs in {key}')

        if manifests[0]['inputs'].get(str(args.reference)) != digest(args.reference):
            raise ValueError('Reference changed since the paired runs')
        truth = structures(args.reference)
        before = structures(args.out / 'baseline' / 'annotation.gff3')
        after = structures(args.out / 'metadata' / 'annotation.gff3')

        summary['changes'] = {
            'unchanged': len(before & after),
            'removed': len(before - after),
            'added': len(after - before),
            'true_removed': len((before - after) & truth),
            'true_added': len((after - before) & truth),
            'false_removed': len((before - after) - truth),
            'false_added': len((after - before) - truth)
        }

    summary_path = args.out / 'summary.json'
    summary_path.write_text(json.dumps(summary, indent=2) + '\n')

    comparison_path = args.out / 'comparison.tsv'
    with comparison_path.open('w') as out:
        writer = csv.DictWriter(out, fieldnames=['variant', 'method', 'level', 'sn', 'pr', 'f1'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--data', type=Path, default=ROOT / 'data' / 'dmel')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--fasta', type=Path)
    p.add_argument('--reference', type=Path)
    p.add_argument('--baseline-binary', type=Path)
    p.add_argument('--metadata-binary', type=Path, default=ROOT / 'src' / 'uniann')
    p.add_argument('--gffread', type=Path, default=ROOT / 'gffread')
    p.add_argument('--min-cds', type=int, default=200)
    p.add_argument('--variants', nargs='*', choices=['baseline', 'metadata'], default=['baseline', 'metadata'])

    args = p.parse_args()
    if args.min_cds < 0:
        p.error('--min-cds must be nonnegative')
    args.data = args.data.resolve()
    args.out = args.out.resolve()

    args.fasta = (args.fasta or args.data / 'NC_004354.4.fa').resolve()
    args.reference = (args.reference or args.data / 'NC_004354.4.p.CDS.uniq.gtf').resolve()
    args.baseline_binary = (args.baseline_binary or ROOT / 'src' / 'uniann').resolve()

    for key in ('metadata_binary', 'gffread'):
        setattr(args, key, getattr(args, key).resolve())

    args.out.mkdir(parents=True, exist_ok=True)
    if args.variants:
        for name in ('summary.json', 'comparison.tsv'):
            (args.out / name).unlink(missing_ok=True)

    inputs = [args.fasta]
    for name in ['out.ps.txt', 'out.gt.txt', 'out.ag.txt', 'out.atg.txt', 'out.stop.txt']:
        inputs.append(args.data / name)

    for variant in args.variants:
        evaluate(args, variant, inputs)

    summarize(args)


if __name__ == '__main__':
    main()
