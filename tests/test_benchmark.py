"""Annotation-comparison fixtures, not measurements of biological accuracy."""
from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from orsearcher.benchmark import annotation_models, evaluate, receptor_family


def feature(kind, start, end, attributes, seqid="chr1", strand="+", phase="."):
    return "\t".join(map(str, (seqid, "fixture", kind, start, end, ".", strand, phase, attributes))) + "\n"


def reference(gene="g1", name="OR1", exons=((1, 900),), seqid="chr1", strand="+", partial=False):
    tid = gene + ".t1"
    text = feature("gene", exons[0][0], exons[-1][1],
                   f"ID={gene};Name={name};partial={str(partial).lower()}", seqid, strand)
    text += feature("mRNA", exons[0][0], exons[-1][1], f"ID={tid};Parent={gene}", seqid, strand)
    count = 0
    coding = []
    for index, (start, end) in enumerate(sorted(exons, reverse=strand == "-"), 1):
        coding.append(feature("CDS", start, end, f"ID={tid}.cds{index};Parent={tid}",
                              seqid, strand, (3 - count % 3) % 3))
        count += end - start + 1
    return text + "".join(coding)


def predictions(folder, records, candidate=False):
    """records: (id, exon pairs, category, optional seqid/strand/phase/independent)."""
    folder.mkdir(parents=True, exist_ok=True)
    stem = "candidate" if candidate else "OR"
    fields = ["candidate_id", "sequence_id", "strand", "category", "confidence", "independent_locus"]
    lines = ["##gff-version 3\n"]
    with (folder / f"{stem}_catalogue.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for record in records:
            cid, exons, category = record[:3]
            opts = record[3] if len(record) > 3 else {}
            seqid, strand = opts.get("seqid", "chr1"), opts.get("strand", "+")
            writer.writerow(dict(candidate_id=cid, sequence_id=seqid, strand=strand,
                                 category=category, confidence="fixture_only",
                                 independent_locus=opts.get("independent", True)))
            if exons:
                lines.append(feature("mRNA", exons[0][0], exons[-1][1], f"ID={cid}", seqid, strand))
            count = 0
            for index, (start, end) in enumerate(sorted(exons, reverse=strand == "-"), 1):
                phase = opts.get("phase", (3 - count % 3) % 3)
                lines.append(feature("CDS", start, end, f"ID={cid}.cds{index};Parent={cid}",
                                     seqid, strand, phase))
                count += end - start + 1
    (folder / f"{stem}_models.gff3").write_text("".join(lines), encoding="utf-8")


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.annotation = self.root / "reference.gff3"
        self.annotation.write_text("##gff-version 3\n" + reference(), encoding="utf-8")

    def evaluate(self, prediction=None, **options):
        args = dict(annotation=self.annotation, prediction=prediction or self.root / "prediction",
                    outdir=self.root / "evaluation")
        args.update(options)
        return evaluate(SimpleNamespace(**args))

    def test_canonical_and_species_prefixed_receptor_names(self):
        for value in ("OR1", "Or67d", "LmigOR1", "CkiaOr123", "Orco", "odorant receptor 1"):
            with self.subTest(value=value):
                self.assertEqual(receptor_family(value), "OR")
        for value in ("GR1", "Gr21a", "LmigGR17", "gustatory receptor 17"):
            with self.subTest(value=value):
                self.assertEqual(receptor_family(value), "GR")
        for value in ("odorant_binding_protein", "ORF", "ORF123", "ordinary_protein"):
            with self.subTest(value=value):
                self.assertEqual(receptor_family(value), "")

    def test_transcript_product_recognizes_generic_gene(self):
        self.annotation.write_text(feature("gene", 1, 900, "ID=g1;Name=LOC123") +
                                   feature("mRNA", 1, 900, "ID=t1;Parent=g1;product=odorant%20receptor%201") +
                                   feature("CDS", 1, 900, "Parent=t1", phase=0), encoding="utf-8")
        genes, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(genes["g1"]["family"], "OR")
        self.assertEqual(len(models), 1)

    def test_gzip_annotation(self):
        compressed = self.root / "reference.gff3.gz"
        with gzip.open(compressed, "wt", encoding="utf-8") as handle:
            handle.write(self.annotation.read_text(encoding="utf-8"))
        genes, _, models, _ = annotation_models(compressed)
        self.assertEqual(set(genes), {"g1"})
        self.assertEqual(len(models), 1)

    def test_direct_gene_parent_for_cds(self):
        self.annotation.write_text(feature("gene", 1, 900, "ID=g1;Name=OR1") +
                                   feature("CDS", 1, 900, "Parent=g1", phase=0), encoding="utf-8")
        genes, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(set(genes), {"g1"})
        self.assertEqual(len(models), 1)

    def test_shared_cds_parents_and_isoforms_count_one_gene(self):
        self.annotation.write_text(feature("gene", 1, 900, "ID=g1;Name=OR1") +
                                   feature("mRNA", 1, 900, "ID=t1;Parent=g1") +
                                   feature("mRNA", 1, 900, "ID=t2;Parent=g1") +
                                   feature("CDS", 1, 900, "Parent=t1,t2", phase=0), encoding="utf-8")
        _, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(len(models), 2)
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        result = self.evaluate()
        self.assertEqual(result["reference_OR_gene_loci"], 1)
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)

    def test_final_results_preferred_to_old_and_unverified_candidates(self):
        sample = self.root / "prediction"
        predictions(sample / "results", [("final", [(1, 900)], "complete")])
        predictions(sample / "nonredundant", [("old", [(2001, 2900)], "partial")])
        predictions(sample, [("unverified", [(2001, 2900)], "partial")], candidate=True)
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["OR_loci_with_complete_splice_endpoint_match"], 1)

    def test_direct_results_folder(self):
        folder = self.root / "prediction" / "results"
        predictions(folder, [("p1", [(1, 900)], "complete")])
        self.assertEqual(self.evaluate(folder)["OR_annotation_loci_recovered"], 1)

    def test_single_batch_sample_is_auto_selected(self):
        predictions(self.root / "prediction" / "samples" / "insect_A" / "results",
                    [("p1", [(1, 900)], "complete")])
        self.assertEqual(self.evaluate()["OR_annotation_loci_recovered"], 1)

    def test_multiple_batch_samples_require_explicit_selection(self):
        for name in ("insect_A", "insect_B"):
            predictions(self.root / "prediction" / "samples" / name / "results",
                        [(name, [(1, 900)], "complete")])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()
        self.assertEqual(self.evaluate(sample="insect_B")["OR_annotation_loci_recovered"], 1)

    def test_incomplete_final_pair_never_silently_uses_older_results(self):
        sample = self.root / "prediction"
        predictions(sample / "results", [("final", [(1, 900)], "complete")])
        (sample / "results" / "OR_models.gff3").unlink()
        predictions(sample / "nonredundant", [("old", [(1, 900)], "complete")])
        with self.assertRaises((ValueError, RuntimeError, FileNotFoundError)):
            self.evaluate()

    def test_explicit_legacy_and_candidate_scopes(self):
        sample = self.root / "prediction"
        predictions(sample / "results", [("final", [(1, 900)], "complete")])
        predictions(sample / "nonredundant", [("old", [(2001, 2900)], "partial")])
        predictions(sample, [("unverified", [(1, 900)], "complete")], candidate=True)
        self.assertEqual(self.evaluate(prediction_scope="legacy", outdir=self.root / "legacy_evaluation")
                         ["OR_annotation_loci_recovered"], 0)
        pool = self.evaluate(prediction_scope="candidates", outdir=self.root / "candidate_evaluation")
        self.assertEqual(pool["OR_annotation_loci_recovered"], 1)
        self.assertEqual(pool["OR_loci_with_complete_splice_endpoint_match"], 0)
        self.assertTrue(pool["records_are_not_gene_counts"])

    def test_zero_recognized_coding_or_is_error_not_perfect_score(self):
        self.annotation.write_text(reference(name="ordinary_protein"), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()

    def test_genes_without_cds_excluded_from_recall_denominator(self):
        with self.annotation.open("a", encoding="utf-8") as handle:
            handle.write(feature("gene", 2001, 2900, "ID=g_no_cds;Name=OR2"))
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        result = self.evaluate()
        self.assertEqual(result["reference_OR_gene_loci"], 1)
        self.assertEqual(result["locus_recall_against_annotation"], 1)

    def test_explicit_gene_family_map(self):
        self.annotation.write_text(reference(name="LOC123"), encoding="utf-8")
        family_map = self.root / "families.tsv"
        family_map.write_text("feature_id\tfamily\tname\ng1\tOR\tcurated_OR1\n", encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        self.assertEqual(self.evaluate(family_map=family_map)["OR_annotation_loci_recovered"], 1)

    def test_family_map_can_match_protein_alias(self):
        self.annotation.write_text(feature("gene", 1, 900, "ID=g1;Name=LOC123") +
                                   feature("mRNA", 1, 900, "ID=t1;Parent=g1") +
                                   feature("CDS", 1, 900, "Parent=t1;protein_id=XP_123.1", phase=0),
                                   encoding="utf-8")
        family_map = self.root / "families.tsv"
        family_map.write_text("feature_id\tfamily\nXP_123.1\tOR\n", encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        self.assertEqual(self.evaluate(family_map=family_map)["OR_annotation_loci_recovered"], 1)

    def test_unused_family_map_entry_is_error(self):
        family_map = self.root / "families.tsv"
        family_map.write_text("feature_id\tfamily\nnot_in_annotation\tOR\n", encoding="utf-8")
        with self.assertRaises((ValueError, RuntimeError)):
            annotation_models(self.annotation, family_map=family_map)

    def test_conflicting_family_map_rows_are_error(self):
        family_map = self.root / "families.tsv"
        family_map.write_text("feature_id\tfamily\ng1\tOR\ng1\tGR\n", encoding="utf-8")
        with self.assertRaises((ValueError, RuntimeError)):
            annotation_models(self.annotation, family_map=family_map)

    def test_conflicting_automatic_gene_and_product_families_are_error(self):
        self.annotation.write_text(feature("gene", 1, 900, "ID=g1;Name=OR1") +
                                   feature("mRNA", 1, 900, "ID=t1;Parent=g1;product=gustatory%20receptor%201") +
                                   feature("CDS", 1, 900, "Parent=t1", phase=0), encoding="utf-8")
        with self.assertRaises((ValueError, RuntimeError)):
            annotation_models(self.annotation)

    def test_explicit_name_pattern_for_nonstandard_ids(self):
        self.annotation.write_text(reference(name="custom_receptor_001"), encoding="utf-8")
        genes, _, _, _ = annotation_models(self.annotation, or_name_patterns=(r"^custom_receptor_\d+$",))
        self.assertEqual(genes["g1"]["family"], "OR")

    def test_missing_prediction_cds_is_error(self):
        predictions(self.root / "prediction", [("p1", [], "complete")])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()

    def test_duplicate_candidate_ids_are_error(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete"),
                                              ("p1", [(1, 900)], "complete")])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()

    def test_disjoint_sequence_ids_require_explicit_mapping(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete", {"seqid": "NC_1.1"})])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()
        seqid_map = self.root / "seqids.tsv"
        seqid_map.write_text("prediction_seqid\tannotation_seqid\nNC_1.1\tchr1\n", encoding="utf-8")
        self.assertEqual(self.evaluate(seqid_map=seqid_map)["OR_annotation_loci_recovered"], 1)

    def test_coordinate_versions_are_not_silently_stripped(self):
        self.annotation.write_text(reference(seqid="NC_1.1"), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete", {"seqid": "NC_1.2"})])
        with self.assertRaises((ValueError, RuntimeError)):
            self.evaluate()

    def test_empty_prediction_is_valid_zero_recall(self):
        predictions(self.root / "prediction", [])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 0)
        self.assertEqual(result["locus_recall_against_annotation"], 0)

    def test_coverage_thresholds_report_separately(self):
        predictions(self.root / "prediction", [("p1", [(1, 600)], "partial")])
        result = self.evaluate()
        self.assertEqual(result["OR_loci_recovered_at_50pct"], 1)
        self.assertEqual(result["OR_loci_recovered_at_75pct"], 0)
        self.assertEqual(result["OR_loci_recovered_at_90pct"], 0)
        self.assertEqual(result["OR_loci_with_complete_splice_endpoint_match"], 0)

    def test_two_predictions_of_one_gene_are_one_recovered_locus(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete"),
                                              ("p2", [(1, 900)], "partial")])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["duplicated_OR_reference_loci"], 1)

    def test_one_fused_prediction_cannot_recover_two_independent_genes(self):
        self.annotation.write_text(reference(exons=((1, 300),)) +
                                   reference(gene="g2", name="OR2", exons=((501, 800),)), encoding="utf-8")
        predictions(self.root / "prediction", [("fusion", [(1, 300), (501, 800)], "partial")])
        result = self.evaluate()
        self.assertEqual(result["reference_OR_gene_loci"], 2)
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["merged_OR_prediction_records"], 1)

    def test_nondependent_damage_alternative_does_not_inflate_primary_recall(self):
        predictions(self.root / "prediction", [("damaged", [(1, 900)], "pseudogene", {"independent": False})])
        self.assertEqual(self.evaluate()["OR_annotation_loci_recovered"], 0)

    def test_gr_overlap_is_not_hidden_by_better_or_hit(self):
        self.annotation.write_text(reference(exons=((1, 900),)) +
                                   reference(gene="g_gr", name="GR1", exons=((601, 900),)), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        result = self.evaluate()
        self.assertEqual(result["predictions_overlapping_GR"], 1)

    def test_strict_structure_match_requires_exact_endpoints(self):
        predictions(self.root / "prediction", [("p1", [(4, 900)], "complete")])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["OR_loci_with_complete_splice_endpoint_match"], 0)

    def test_strict_structure_match_requires_cds_phase(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete", {"phase": 1})])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["OR_loci_with_complete_splice_endpoint_match"], 0)

    def test_strict_structure_match_rejects_partial_reference(self):
        self.annotation.write_text(reference(partial=True), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 1)
        self.assertEqual(result["OR_loci_with_complete_splice_endpoint_match"], 0)
        self.assertEqual(result["OR_loci_with_coordinate_phase_match"], 1)
        self.assertEqual(result["OR_loci_with_complete_prediction_coordinate_phase_match"], 1)

    def test_minus_strand_phase_and_splice_match(self):
        exons = ((1, 298), (401, 1002))
        self.annotation.write_text(reference(exons=exons, strand="-"), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", list(exons), "complete", {"strand": "-"})])
        self.assertEqual(self.evaluate()["OR_loci_with_complete_splice_endpoint_match"], 1)

    def test_opposite_strand_overlap_does_not_recover_gene(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete", {"strand": "-"})])
        self.assertEqual(self.evaluate()["OR_annotation_loci_recovered"], 0)

    def test_matching_reassigns_ambiguous_prediction_to_recover_both_loci(self):
        self.annotation.write_text(reference(exons=((1, 1000),)) +
                                   reference(gene="g2", name="OR2", exons=((1001, 2000),)), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 2000)], "partial"),
                                              ("p2", [(1, 500)], "partial")])
        result = self.evaluate()
        self.assertEqual(result["OR_annotation_loci_recovered"], 2)
        self.assertEqual(result["merged_OR_prediction_records"], 1)

    def test_overlapping_cds_uses_union_coverage_not_double_counting(self):
        self.annotation.write_text(reference(exons=((1, 1000),)), encoding="utf-8")
        predictions(self.root / "prediction", [("p1", [(1, 400), (301, 600)], "partial")])
        result = self.evaluate()
        self.assertEqual(result["OR_loci_recovered_at_50pct"], 1)
        self.assertEqual(result["OR_loci_recovered_at_75pct"], 0)
        self.assertEqual(result["prediction_records_with_overlapping_CDS"], 1)
        with (self.root / "evaluation" / "all_annotation_overlaps.tsv").open(encoding="utf-8") as handle:
            record = next(csv.DictReader(handle, delimiter="\t"))
        self.assertAlmostEqual(float(record["reference_CDS_covered"]), 0.6)

    def test_one_map_alias_cannot_silently_label_two_genes(self):
        text = (reference(name="LOC123").replace("Name=LOC123", "Name=LOC123;Alias=shared") +
                reference(gene="g2", name="LOC456", exons=((2001, 2900),))
                .replace("Name=LOC456", "Name=LOC456;Alias=shared"))
        self.annotation.write_text(text, encoding="utf-8")
        family_map = self.root / "families.tsv"
        family_map.write_text("feature_id\tfamily\nshared\tOR\n", encoding="utf-8")
        with self.assertRaises((ValueError, RuntimeError)):
            annotation_models(self.annotation, family_map=family_map)

    def test_evaluation_refuses_to_overwrite_completed_output(self):
        predictions(self.root / "prediction", [("p1", [(1, 900)], "complete")])
        self.evaluate()
        result_path = self.root / "evaluation" / "benchmark_summary.json"
        original = result_path.read_bytes()
        with self.assertRaises(FileExistsError):
            self.evaluate()
        self.assertEqual(result_path.read_bytes(), original)

    def test_embedded_fasta_is_not_parsed_as_gff(self):
        with self.annotation.open("a", encoding="utf-8") as handle:
            handle.write("##FASTA\n>chr1\nATGCATGC\n")
        genes, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(set(genes), {"g1"})
        self.assertEqual(len(models), 1)

    def test_id_only_receptor_symbol_is_recognized(self):
        self.annotation.write_text(feature('gene', 1, 900, 'ID=LmigOR1') +
                                   feature('CDS', 1, 900, 'Parent=LmigOR1', phase=0), encoding='utf-8')
        genes, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(genes['LmigOR1']['family'], 'OR')
        self.assertEqual(len(models), 1)

    def test_transcript_uncertain_endpoint_is_not_complete_truth(self):
        text = reference().replace('ID=g1.t1;Parent=g1', 'ID=g1.t1;Parent=g1;start_range=.,1')
        self.annotation.write_text(text, encoding='utf-8')
        predictions(self.root / 'prediction', [('p1', [(1, 900)], 'complete')])
        self.assertEqual(self.evaluate()['OR_loci_with_complete_splice_endpoint_match'], 0)

    def test_gene_pseudo_flag_is_not_intact_complete_truth(self):
        text = reference().replace('Name=OR1', 'Name=OR1;pseudo=true')
        self.annotation.write_text(text, encoding='utf-8')
        predictions(self.root / 'prediction', [('p1', [(1, 900)], 'complete')])
        result = self.evaluate()
        self.assertEqual(result['OR_annotation_loci_recovered'], 1)
        self.assertEqual(result['OR_loci_with_complete_splice_endpoint_match'], 0)

    def test_prediction_on_nonreceptor_annotation_contig_is_unmatched_not_invalid(self):
        text = reference() + reference(gene='other', name='actin', seqid='chr2')
        self.annotation.write_text(text, encoding='utf-8')
        predictions(self.root / 'prediction', [('p1', [(1, 900)], 'complete', {'seqid': 'chr2'})])
        result = self.evaluate()
        self.assertEqual(result['OR_annotation_loci_recovered'], 0)
        self.assertEqual(result['annotation_unmatched_predictions'], 1)

    def test_catalogue_bounds_disagreeing_with_gff_are_error(self):
        folder = self.root / 'prediction'
        predictions(folder, [('p1', [(1, 900)], 'complete')])
        path = folder / 'OR_catalogue.tsv'
        lines = path.read_text(encoding='utf-8').splitlines()
        path.write_text(lines[0] + '\tstart\tend\n' + lines[1] + '\t4\t900\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'bounds mismatch'):
            self.evaluate()

    def test_benchmark_cli_options_are_available_without_prediction_dependencies(self):
        from orsearcher.cli import parser
        args = parser().parse_args(['benchmark', '--annotation', 'gold.gff3', '--prediction', 'run',
                                    '--outdir', 'score', '--sample', 'species_A', '--prediction-scope', 'final',
                                    '--family-map', 'families.tsv', '--or-name-pattern', '^custom_OR',
                                    '--seqid-map', 'contigs.tsv', '--min-reference-coverage', '0.75'])
        self.assertEqual(args.sample, 'species_A')
        self.assertEqual(args.min_reference_coverage, .75)
        self.assertEqual(args.or_name_pattern, ['^custom_OR'])

    def test_unrelated_trans_spliced_multipart_gene_does_not_block_or_inventory(self):
        text = reference()
        text += feature('gene', 1001, 1300, 'ID=complex;Name=mod;part=1;exception=trans-splicing', strand='+')
        text += feature('gene', 1401, 1700, 'ID=complex;Name=mod;part=2;exception=trans-splicing', strand='-')
        text += feature('mRNA', 1401, 1700, 'ID=complex.t;Parent=complex', strand='-')
        text += feature('CDS', 1401, 1700, 'Parent=complex.t', strand='-', phase=0)
        self.annotation.write_text(text, encoding='utf-8')
        genes, _, models, _ = annotation_models(self.annotation)
        self.assertEqual(set(genes), {'g1'})
        self.assertEqual(len(models), 1)

    def test_selected_receptor_multipart_on_opposite_strands_is_explicit_error(self):
        text = feature('gene', 1, 300, 'ID=complex;Name=OR1;part=1;exception=trans-splicing', strand='+')
        text += feature('gene', 401, 700, 'ID=complex;Name=OR1;part=2;exception=trans-splicing', strand='-')
        text += feature('CDS', 1, 300, 'Parent=complex', phase=0)
        self.annotation.write_text(text, encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'unsupported multi-contig/strand'):
            annotation_models(self.annotation)


if __name__ == "__main__":
    unittest.main()
