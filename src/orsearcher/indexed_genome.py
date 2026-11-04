"""Bounded random access for plain indexed FASTA, samtools for compressed input."""
from collections import OrderedDict
from pathlib import Path
import subprocess


class IndexedGenome:
    def __init__(self, path, cache_bytes=8 * 1024 * 1024):
        self.path = Path(path)
        self.handle = None
        self.index = {}
        self.cache = OrderedDict()
        self.cache_bytes = cache_bytes
        self.used = 0

    def __enter__(self):
        index = Path(str(self.path) + '.fai')
        with self.path.open('rb') as handle:
            compressed = handle.read(2) == b'\x1f\x8b'
        if index.exists() and not compressed:
            for line in index.read_text().splitlines():
                fields = line.split('\t')
                name = fields[0]
                length, offset, bases, width = map(int, fields[1:5])
                if name in self.index or length < 0 or offset < 0 or bases <= 0 or width < bases:
                    raise ValueError('Invalid FASTA index: ' + name)
                self.index[name] = (length, offset, bases, width)
            self.handle = self.path.open('rb')
        return self

    def __exit__(self, *args):
        if self.handle is not None:
            self.handle.close()
        self.cache.clear()

    def fetch(self, contig, start, end):
        if start < 1 or end < start:
            raise ValueError('Invalid genome interval')
        key = (contig, start, end)
        if key in self.cache:
            self.cache.move_to_end(key)
            return self.cache[key]
        if self.handle is None:
            result = subprocess.run(['samtools', 'faidx', str(self.path), f'{contig}:{start}-{end}'],
                                    capture_output=True, text=True, check=True)
            sequence = ''.join(result.stdout.splitlines()[1:]).upper()
        else:
            length, offset, bases, width = self.index[contig]
            if end > length:
                raise ValueError('Genome interval exceeds contig: ' + contig)
            first = offset + (start - 1) // bases * width + (start - 1) % bases
            last = offset + (end - 1) // bases * width + (end - 1) % bases
            self.handle.seek(first)
            raw = self.handle.read(last - first + 1)
            sequence = raw.replace(b'\n', b'').replace(b'\r', b'').decode('ascii').upper()
        if len(sequence) != end - start + 1:
            raise ValueError('Genome extraction length mismatch: ' + contig)
        if len(sequence) <= self.cache_bytes:
            while self.cache and self.used + len(sequence) > self.cache_bytes:
                _, old = self.cache.popitem(last=False)
                self.used -= len(old)
            self.cache[key] = sequence
            self.used += len(sequence)
        return sequence
