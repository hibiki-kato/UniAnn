"""Check exact CDS matching, the wrapper filter, and resource measurements."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('compare_viterbi', ROOT / 'benchmarks/compare_viterbi.py')
benchmark = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(benchmark)


def cds(start, end, phase, attributes):
    return f'chr1\ttest\tCDS\t{start}\t{end}\t.\t+\t{phase}\t{attributes}\n'


class BenchmarkTests(unittest.TestCase):
    def test_gff3_and_gtf_match_with_shared_and_duplicate_transcripts(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / 'reference.gtf'
            query = Path(directory) / 'query.gff3'
            reference.write_text(cds(1, 9, 0, 'transcript_id "a";') +
                                 cds(20, 25, 0, 'transcript_id "a";'))
            query.write_text(cds(20, 25, 0, 'Parent=b,c') + cds(1, 9, 0, 'Parent=b,c'))
            result = benchmark.exact_metrics(reference, query)
            self.assertEqual((result['reference'], result['predictions'], result['tp']), (1, 1, 1))
            self.assertEqual(result['f1'], 100)
            query.write_text(cds(20, 25, 1, 'Parent=b') + cds(1, 9, 0, 'Parent=b'))
            result = benchmark.exact_metrics(reference, query)
            self.assertEqual((result['tp'], result['fp'], result['fn']), (0, 1, 1))

    def test_missing_parent_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            query = Path(directory) / 'query.gff3'
            for attributes in ('ID=cds1', 'Parent='):
                query.write_text(cds(1, 9, 0, attributes))
                with self.assertRaises(ValueError):
                    benchmark.structures(query)

    def test_filter_threshold_and_exon_count(self):
        # gffread output has CDS-only exons for the decoder's raw predictions.
        rows = [
            'exonCount=1;exons=1-200;CDS=1:200;CDSphase=0',
            'exonCount=2;exons=1-100,201-301;CDS=1:301;CDSphase=0',
            'exonCount=3;exons=1-3,5-7,9-11;CDS=1:11;CDSphase=0',
        ]
        lines = [f'chr1\ttest\ttranscript\t1\t301\t.\t+\t.\tID=t{i};{row}\n'
                 for i, row in enumerate(rows)]
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            def fake_execute(command, stdout, stderr):
                if '--tlf' in command:
                    stdout.write_text(''.join(lines))
                else:
                    stdout.write_text(Path(command[1]).read_text())
            with patch.object(benchmark, 'execute', side_effect=fake_execute):
                annotation = benchmark.postprocess(directory / 'raw', directory, Path('gffread'), 200)
            self.assertEqual(annotation.read_text(), ''.join(lines[1:]))

    def test_resource_measurements(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / 'time.log'
            for elapsed, seconds in [('0:41.55', 41.55), ('1:02:03.5', 3723.5)]:
                log.write_text(f'Elapsed (wall clock) time (h:mm:ss or m:ss): {elapsed}\n'
                               'Maximum resident set size (kbytes): 1406960\n')
                self.assertEqual(benchmark.parse_time_v(log),
                                 {'wall_seconds': seconds, 'max_rss_kib': 1406960})
            log.write_text('')
            with self.assertRaises(ValueError):
                benchmark.parse_time_v(log)


if __name__ == '__main__':
    unittest.main()
