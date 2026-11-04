from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from orsearcher.models import (
    attach_reference_lengths,
    classify_candidate,
    deduplicate_candidates,
    extract_sequences,
    parse_miniprot_gff,
)
from orsearcher.sharding import Segment


class ModelTests(unittest.TestCase):
    def make_mapping(self) -> dict[str, Segment]:
        return {
            "ors0000000001": Segment(
                "ors0000000001",
                "chr1",
                "chr1 chromosome",
                1000,
                1300,
                5000,
                "shard_00001.fa",
                True,
            )
        }

    def test_parse_remap_extract_and_classify(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hit.gff3"
            path.write_text(
                "##gff-version 3\n"
                "ors0000000001\tminiprot\tmRNA\t10\t108\t.\t+\t.\t"
                "ID=MP000001;Rank=1;Identity=0.45;Positive=0.60;"
                "Frameshift=0;StopCodon=0;Target=ref001 1 32\n"
                "ors0000000001\tminiprot\tCDS\t10\t105\t.\t+\t0\tParent=MP000001\n"
                "ors0000000001\tminiprot\tstop_codon\t106\t108\t.\t+\t0\tParent=MP000001\n",
                encoding="utf-8",
            )
            candidates = parse_miniprot_gff(path, self.make_mapping())
            attach_reference_lengths(candidates, {"ref001": 32})
            candidate = candidates[0]
            self.assertEqual((candidate.original_start, candidate.original_end), (1010, 1108))
            sequence = "N" * 9 + "ATG" + "GCT" * 30 + "TAA" + "N" * (300 - 105)
            extract_sequences(candidate, sequence)
            classify_candidate(candidate, 0.15, 0.25, 20, 100)
            self.assertEqual(candidate.status, "functional_candidate")
            self.assertTrue(candidate.has_start)
            self.assertTrue(candidate.has_stop)
            self.assertEqual(len(candidate.protein_sequence), 31)

    def test_deduplicate_retains_best_model_and_reference_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hits.gff3"
            path.write_text(
                "ors0000000001\tminiprot\tmRNA\t10\t108\t.\t+\t.\t"
                "ID=A;Rank=2;Identity=0.30;Positive=0.40;Frameshift=0;StopCodon=0;Target=refA 1 30\n"
                "ors0000000001\tminiprot\tCDS\t10\t108\t.\t+\t0\tParent=A\n"
                "ors0000000001\tminiprot\tmRNA\t12\t109\t.\t+\t.\t"
                "ID=B;Rank=1;Identity=0.50;Positive=0.60;Frameshift=0;StopCodon=0;Target=refB 1 32\n"
                "ors0000000001\tminiprot\tCDS\t12\t109\t.\t+\t0\tParent=B\n",
                encoding="utf-8",
            )
            candidates = parse_miniprot_gff(path, self.make_mapping())
            attach_reference_lengths(candidates, {"refA": 40, "refB": 32})
            kept = deduplicate_candidates(candidates)
            self.assertEqual(len(kept), 1)
            self.assertEqual(kept[0].reference_id, "refB")
            self.assertEqual(kept[0].reference_hits, {"refA", "refB"})

    def test_deduplicate_does_not_bridge_tandem_genes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "tandem.gff3"
            path.write_text(
                "ors0000000001\tminiprot\tmRNA\t10\t108\t.\t+\t.\t"
                "ID=A;Rank=1;Identity=0.95;Positive=0.98;Target=refA 1 33\n"
                "ors0000000001\tminiprot\tCDS\t10\t108\t.\t+\t0\tParent=A\n"
                "ors0000000001\tminiprot\tmRNA\t210\t308\t.\t+\t.\t"
                "ID=B;Rank=1;Identity=0.96;Positive=0.99;Target=refB 1 33\n"
                "ors0000000001\tminiprot\tCDS\t210\t308\t.\t+\t0\tParent=B\n"
                "ors0000000001\tminiprot\tmRNA\t10\t308\t.\t+\t.\t"
                "ID=BRIDGE;Rank=3;Identity=0.30;Positive=0.45;Target=refBridge 1 33\n"
                "ors0000000001\tminiprot\tCDS\t10\t308\t.\t+\t0\tParent=BRIDGE\n",
                encoding="utf-8",
            )
            candidates = parse_miniprot_gff(path, self.make_mapping())
            attach_reference_lengths(
                candidates, {"refA": 33, "refB": 33, "refBridge": 33}
            )
            kept = deduplicate_candidates(candidates)
            self.assertEqual(len(kept), 2)
            self.assertEqual({item.reference_id for item in kept}, {"refA", "refB"})
            self.assertIn(
                "refBridge",
                set().union(*(item.reference_hits for item in kept)),
            )

    def test_disruption_is_pseudogene(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "hit.gff3"
            path.write_text(
                "ors0000000001\tminiprot\tmRNA\t10\t108\t.\t+\t.\t"
                "ID=P;Rank=1;Identity=0.40;Positive=0.50;Frameshift=1;StopCodon=0;Target=ref 1 30\n"
                "ors0000000001\tminiprot\tCDS\t10\t108\t.\t+\t0\tParent=P\n",
                encoding="utf-8",
            )
            candidate = parse_miniprot_gff(path, self.make_mapping())[0]
            attach_reference_lengths([candidate], {"ref": 30})
            extract_sequences(candidate, "N" * 9 + "ATG" + "GCT" * 31 + "TAA" + "N" * 195)
            classify_candidate(candidate, 0.15, 0.25, 20, 100)
            self.assertEqual(candidate.status, "pseudogene_candidate")

    def test_multiexon_phase_does_not_delete_split_codon(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "spliced.gff3"
            path.write_text(
                "ors0000000001\tminiprot\tmRNA\t10\t25\t.\t+\t.\t"
                "ID=S;Rank=0;Identity=0.8;Positive=0.9;Target=ref 1 2\n"
                "ors0000000001\tminiprot\tCDS\t10\t14\t.\t+\t0\tParent=S\n"
                "ors0000000001\tminiprot\tCDS\t25\t25\t.\t+\t1\tParent=S\n",
                encoding="utf-8",
            )
            candidate = parse_miniprot_gff(path, self.make_mapping())[0]
            # Exon 1 is ATGGC; exon 2 supplies T, together encoding ATG GCT.
            sequence = "N" * 9 + "ATGGC" + "N" * 10 + "T" + "N" * 275
            extract_sequences(candidate, sequence)
            self.assertEqual(candidate.cds_sequence, "ATGGCT")
            self.assertEqual(candidate.protein_sequence, "MA")


if __name__ == "__main__":
    unittest.main()

