"""Exercise the public flags, valid small inputs, and rejected input errors."""
from pathlib import Path
import json
import subprocess
import tempfile
import unittest

from fixtures import cases, write_inputs

BINARY = Path(__file__).resolve().parents[1] / 'src' / 'uniann'


class MetadataCliTests(unittest.TestCase):
    def test_metadata_mode(self):
        for name, sequence, random_scores in cases():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                files = write_inputs(directory, sequence, random_scores)
                command = [str(BINARY), *map(str, files), '--metadata-state-viterbi']
                first = subprocess.run(command, capture_output=True, check=True)
                second = subprocess.run(command + ['--no-dp-dump'], capture_output=True, check=True)
                self.assertTrue(first.stdout.startswith(b'##gff-version 3\n'))
                self.assertEqual(first.stdout, second.stdout)
                self.assertEqual(first.stderr, b'')
                self.assertEqual(second.stderr, b'')

    def test_unchanged_paths_keep_upstream_gff_scores(self):
        expected = json.loads((BINARY.parents[1] / 'tests' / 'legacy_expected.json').read_text())
        for name, sequence, random_scores in cases():
            if name not in ('single_exon', 'early_start', 'intron_40', 'intron_41'):
                continue
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                files = write_inputs(directory, sequence, random_scores)
                result = subprocess.run([str(BINARY), *map(str, files), '--metadata-state-viterbi'],
                                        capture_output=True, check=True)
                self.assertEqual(result.stdout.decode(), expected[name]['gff'])

    def test_legacy_dump_is_optional(self):
        with tempfile.TemporaryDirectory() as directory:
            _, seq, random_scores = next(cases())
            files = write_inputs(directory, seq, random_scores)
            command = [str(BINARY), *map(str, files)]
            normal = subprocess.run(command, capture_output=True, check=True)
            quiet = subprocess.run(command + ['--no-dp-dump'], capture_output=True, check=True)
            self.assertEqual(normal.stdout, quiet.stdout)
            self.assertIn(b'\tdp\t', normal.stderr)
            self.assertEqual(quiet.stderr, b'')

    def test_invalid_arguments_and_empty_sequence(self):
        missing = subprocess.run([str(BINARY)], capture_output=True)
        self.assertEqual(missing.returncode, 1)
        self.assertIn(b'atg.txt stop.txt', missing.stderr)
        with tempfile.TemporaryDirectory() as directory:
            files = write_inputs(directory, '')
            result = subprocess.run([str(BINARY), *map(str, files), '--metadata-state-viterbi'], capture_output=True)
            self.assertEqual(result.returncode, 1)
            unknown = subprocess.run([str(BINARY), *map(str, files), '--unknown'], capture_output=True)
            self.assertEqual(unknown.returncode, 1)
            self.assertIn(b'Unknown option', unknown.stderr)

    def test_invalid_input_positions(self):
        for filename, bad_line in [('emissions.txt', '-1 0 0 0 0 0\n'),
                                   ('emissions.txt', '0 0 1\n'), ('gt.txt', '9 1\n')]:
            with self.subTest(filename=filename, line=bad_line), tempfile.TemporaryDirectory() as directory:
                files = write_inputs(directory, 'AAA')
                (Path(directory) / filename).write_text(bad_line)
                result = subprocess.run([str(BINARY), *map(str, files), '--metadata-state-viterbi'], capture_output=True)
                self.assertEqual(result.returncode, 1)
                self.assertIn(b'Invalid', result.stderr)


if __name__ == '__main__':
    unittest.main()
