import random
import shutil
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from orsearcher.indexed_genome import IndexedGenome


class IndexedTests(unittest.TestCase):
    def test_wrapping_crlf_and_unwrapped_boundaries(self):
        seq = 'AcGtNRy' * 33
        for width in (1, 7, 80, len(seq)):
            for newline in (b'\n', b'\r\n'):
                with TemporaryDirectory() as tmp:
                    p = Path(tmp) / 'genome.fa'
                    header = b'>chr1 description' + newline
                    p.write_bytes(header + newline.join(seq[i:i+width].encode() for i in range(0,len(seq),width)) + newline)
                    Path(str(p)+'.fai').write_text(f'chr1\t{len(seq)}\t{len(header)}\t{width}\t{width+len(newline)}\n')
                    with IndexedGenome(p, cache_bytes=50) as reader:
                        for start,end in [(1,1),(1,len(seq)),(len(seq),len(seq)),(6,19),(7,7),(8,21)]:
                            self.assertEqual(reader.fetch('chr1',start,end),seq[start-1:end].upper())
                        self.assertLessEqual(reader.used,50)
                        with self.assertRaises(ValueError): reader.fetch('chr1',1,len(seq)+1)
                        with self.assertRaises(KeyError): reader.fetch('missing',1,2)
                        with self.assertRaises(ValueError): reader.fetch('chr1',0,2)

    @unittest.skipUnless(shutil.which('samtools'),'requires native samtools')
    def test_random_intervals_match_samtools(self):
        rng=random.Random(31)
        with TemporaryDirectory() as tmp:
            p=Path(tmp)/'genome.fa'
            p.write_text('>a\n'+ '\n'.join(''.join(rng.choices('ACGTN',k=71)) for _ in range(20))+'\n>b\nACGT\n')
            subprocess.run(['samtools','faidx',str(p)],check=True)
            with IndexedGenome(p) as reader:
                for _ in range(40):
                    start=rng.randint(1,1400);end=rng.randint(start,1420)
                    result=subprocess.run(['samtools','faidx',str(p),f'a:{start}-{end}'],capture_output=True,text=True,check=True)
                    self.assertEqual(reader.fetch('a',start,end),''.join(result.stdout.splitlines()[1:]))
                self.assertEqual(reader.fetch('b',1,4),'ACGT')
