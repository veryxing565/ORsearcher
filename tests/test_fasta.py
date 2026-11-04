from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orsearcher.fasta import detect_sequence_type, longest_orf, read_fasta, reverse_complement, translate


class FastaTests(unittest.TestCase):
    def test_translation_and_reverse_complement(self) -> None:
        self.assertEqual(translate("ATGGCTTAA"), "MA*")
        self.assertEqual(reverse_complement("ATGCRY"), "RYGCAT")

    def test_detect_and_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            dna = root / "or.fna"
            # 120 aa ORF, followed by a stop.
            dna.write_text(">or one\nATG" + "GCT" * 119 + "TAA\n", encoding="utf-8")
            self.assertEqual(detect_sequence_type(("ATG" + "GCT" * 119 + "TAA")), "dna")
            records = list(read_fasta(dna))
            self.assertEqual(len(records),1)
            self.assertEqual(translate(records[0].sequence),'M'+'A'*119+'*')

    def test_longest_orf_can_use_reverse_strand(self) -> None:
        coding = "ATG" + "GCT" * 109 + "TAA"
        protein, strand = longest_orf(reverse_complement(coding))
        self.assertGreaterEqual(len(protein), 110)
        self.assertTrue(strand.startswith("-"))


if __name__ == "__main__":
    unittest.main()

