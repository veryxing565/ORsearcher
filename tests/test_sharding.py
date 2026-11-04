from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orsearcher.fasta import read_fasta
from orsearcher.sharding import create_shards


class ShardingTests(unittest.TestCase):
    def test_long_record_overlap_covers_boundaries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "genome.fa"
            genome.write_text(">chromosome\n" + "A" * 250 + "\n", encoding="utf-8")
            result = create_shards(genome, root / "shards", 100, 30, 250, resume=False)
            coordinates = [(segment.original_start, segment.original_end) for segment in result.segments]
            self.assertEqual(coordinates, [(0, 100), (70, 170), (140, 240), (210, 250)])
            # Any interval no longer than overlap is fully contained in a window.
            for start in range(0, 221):
                end = start + 30
                self.assertTrue(any(left <= start and end <= right for left, right in coordinates))

    def test_small_contigs_are_packed_not_one_file_each(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "fragmented.fa"
            genome.write_text(
                "".join(f">contig_{index}\n{'ACGT' * 5}\n" for index in range(100)),
                encoding="utf-8",
            )
            result = create_shards(genome, root / "shards", 100, 20, 500, resume=False)
            self.assertEqual(result.sequence_count, 100)
            self.assertEqual(len(result.segments), 100)
            self.assertEqual(len(result.shards), 4)
            self.assertEqual(sum(1 for path in result.shards for _ in read_fasta(path)), 100)

    def test_resume_reuses_valid_shards(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            genome = root / "genome.fa"
            genome.write_text(">c\n" + "ACGT" * 50 + "\n", encoding="utf-8")
            first = create_shards(genome, root / "shards", 100, 20, 500, resume=True)
            mtimes = [path.stat().st_mtime_ns for path in first.shards]
            second = create_shards(genome, root / "shards", 100, 20, 500, resume=True)
            self.assertEqual(mtimes, [path.stat().st_mtime_ns for path in second.shards])


if __name__ == "__main__":
    unittest.main()


