"""Check wrapper flag forwarding without invoking external gene-model tools."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


INPUT_ARGS = ['sequence.fa', 'out.ps.txt', 'out.gt.txt', 'out.ag.txt', 'out.atg.txt', 'out.stop.txt']


def decoder_args(wrapper_flags):
    """Run the wrapper with stub tools and return the arguments the decoder saw."""
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        shutil.copyfile(ROOT / 'scripts' / 'uniann.sh', work / 'uniann.sh')
        (work / 'uniann').write_text('#!/bin/sh\nprintf "%s\\n" "$@" > decoder.args\n')
        (work / 'gffread').write_text('#!/bin/sh\ncat\n')
        for name in ('uniann', 'gffread'):
            (work / name).chmod(0o755)
        (work / 'sequence.fa').write_text('>fixture\nGT\n')
        (work / 'psauron.csv').write_text('unused: emission input is already prepared\n')
        (work / 'out.ps.txt').write_text('0 0 0 0 0 0\n')
        (work / 'sites.tsv').write_text('fixture\t1\t+\tdonor\tGT\t0.99\n')
        env = dict(os.environ, PATH=str(work) + os.pathsep + os.environ['PATH'])
        subprocess.run(['bash', 'uniann.sh', *wrapper_flags,
                        '-f', 'sequence.fa', '-p', 'psauron.csv', '-s', 'sites.tsv'],
                       cwd=work, env=env, capture_output=True, check=True, start_new_session=True)
        return (work / 'decoder.args').read_text().splitlines()


class WrapperTests(unittest.TestCase):
    def test_noviterbi_skips_dump_and_keeps_following_flag(self):
        self.assertEqual(decoder_args(['-n', '-m', '2']), INPUT_ARGS + ['--no-dp-dump'])

    def test_default_keeps_dump(self):
        self.assertEqual(decoder_args([]), INPUT_ARGS)


if __name__ == '__main__':
    unittest.main()
