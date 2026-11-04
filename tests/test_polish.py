import copy
import csv
import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from orsearcher import polish
from orsearcher.fasta import FastaRecord


def model(name='a', start=100, full=True, reference='core'):
    m = dict(id=name, cid=name, contig='chr1', strand='+', exons=[[start, start + 1199]],
             cds='ATG' + 'GCT' * 398 + 'TAA', protein='M' + 'A' * 398, length=399,
             has_start=full, has_stop=full, frameshifts=0, stops=[], ambiguous=0, remainder=0,
             canonical=True, bridge=False, fusion_pairs=[], or_bits=700, gr_bits=0, family='OR',
             repeated_or_domains=False, repeated_reference_regions=False,
             or_profile_domains=[dict(domain_evalue=1e-50, coverage=.95)],
             homology_hits=[dict(ref=reference, evalue=1e-150, bits=700, identity=90,
                                 qs=1, qe=399, qlen=399, ss=1, se=399, slen=399, curated=True)],
             _group_ids={name}, _tools={'miniprot'}, _baseline_selected=True)
    m['_assessment'] = polish.assessment(m, {'core'})
    return m


class PolicyTests(unittest.TestCase):
    def test_template_beyond_eight_alternatives_is_used_for_repair(self):
        cluster=[model(str(i),full=False,reference='partial_ref') for i in range(9)]
        cluster.append(model('last',full=False))
        jobs,deferred,unavailable=polish.repair_plan([cluster],{}, {'core':FastaRecord('core','','M'+'A'*398)}, {'core'},SimpleNamespace())
        self.assertEqual(jobs[0]['references'],['core'])
        self.assertFalse(unavailable)

    def test_complete_blockers_distinguish_template_from_missing_stop(self):
        m=model(reference='partial_ref')
        self.assertEqual(m['_assessment']['completeness_blockers'],['no_qualified_full_length_template_support'])
        m['has_stop']=False
        self.assertIn('missing_stop',polish.assessment(m,{'core'})['completeness_blockers'])

    def test_single_tool_complete_can_be_usable(self):
        self.assertTrue(model()['_assessment']['usable_complete'])

    def test_partial_reference_does_not_prove_full_length(self):
        self.assertFalse(model(reference='partial_ref')['_assessment']['usable_complete'])

    def test_damage_and_noncanonical_are_not_usable(self):
        for field, value in [('frameshifts', 1), ('stops', [1]), ('canonical', False), ('ambiguous', 1)]:
            m = model()
            m[field] = value
            self.assertFalse(polish.assessment(m, {'core'})['usable_complete'])

    def test_gr_competition_not_accepted_even_with_cached_positive_label(self):
        m = model()
        m.update(preliminary_identity={'status': 'OR_supported'}, gr_homology_hits=[dict(
            ref='gr1', evalue=1e-150, bits=750, identity=95, qs=1, qe=399, qlen=399,
            ss=1, se=399, slen=399)], gr_profile_domains=[[dict(domain_evalue=1e-50, coverage=.95)]])
        self.assertFalse(polish.assessment(m, {'core'})['identity_ok'])

    def test_background_conflict_is_not_complete(self):
        m = model()
        m['background_hits'] = [dict(ref='other', evalue=1e-150, bits=800, identity=95,
                                    qs=1, qe=399, qlen=399, ss=1, se=399, slen=399)]
        self.assertFalse(polish.assessment(m, {'core'})['usable_complete'])

    def test_identical_peptides_different_loci_stay_separate(self):
        self.assertEqual(len(polish.consolidate([model('a'), model('b', 3000)])[0]), 2)

    def test_contained_frame_compatible_fragment_alias(self):
        a, b = model(), model('partial', full=False)
        b.update(exons=[[400, 999]], cds='GCT' * 200, protein='A' * 200, length=200)
        b['_assessment'] = dict(a['_assessment'], usable_complete=False, credible_partial=True)
        self.assertTrue(polish.same_locus(a, b))
        self.assertEqual(len(polish.consolidate([a, b])[0]), 1)

    def test_wrong_frame_is_not_automatic_alias(self):
        a, b = model(), model('partial', full=False)
        b.update(exons=[[401, 1000]], cds='GCT' * 200, protein='A' * 200, length=200)
        self.assertFalse(polish.same_locus(a, b))

    def test_shared_splice_alternatives_are_one_locus(self):
        a = model()
        a['exons'] = [[100, 399], [1000, 1899]]
        b = copy.deepcopy(a)
        b.update(id='b', exons=[[103, 399], [1000, 1899]])
        self.assertTrue(polish.same_locus(a, b))
        self.assertEqual(len(polish.consolidate([a, b])[0]), 1)

    def test_bridge_does_not_merge_tandem_genes(self):
        a, b, bridge = model('left', 100), model('right', 2200), model('bridge', 100)
        bridge['exons'] = [[100, 699], [2200, 2799]]
        bridge['_assessment'] = dict(bridge['_assessment'], bits=100)
        clusters, assignments, _ = polish.consolidate([a, b, bridge])
        self.assertEqual(len(clusters), 2)
        self.assertTrue(any(r['disposition'] == 'possible_two_gene_fusion' for r in assignments))

    def test_opposite_strands_stay_separate(self):
        a, b = model('a'), model('b')
        b['strand'] = '-'
        self.assertEqual(len(polish.consolidate([a, b])[0]), 2)

    def test_small_coding_overlap_held_not_erased(self):
        a, b = model('a', 100), model('b', 1260)
        clusters, _, conflicts = polish.consolidate([a, b])
        self.assertEqual(len(clusters), 1)
        self.assertEqual(conflicts[0]['overlap_nt'], 40)

    def test_shared_span_without_coding_overlap_not_merged(self):
        a, b = model('a'), model('b')
        a['exons'] = [[100, 399], [2000, 2899]]
        b['exons'] = [[500, 1699]]
        self.assertFalse(polish.same_locus(a, b))
        self.assertEqual(len(polish.consolidate([a, b])[0]), 2)

    def test_minus_frame_check(self):
        a, b = model('a'), model('b')
        a['strand'] = b['strand'] = '-'
        b['exons'] = [[400, 1299]]
        self.assertTrue(polish.same_frame(a, b))
        b['exons'] = [[400, 1298]]
        self.assertFalse(polish.same_frame(a, b))


class InterfaceTests(unittest.TestCase):
    def test_active_complete_ids_respect_holdout(self):
        refs = {'kept': FastaRecord('kept', '', 'M' * 300)}
        ids, _ = polish.complete_reference_ids(SimpleNamespace(complete_reference_ids={'held', 'kept'}, complete_reference='/unread/path'), refs)
        self.assertEqual(ids, {'kept'})

    def test_no_complete_metadata_does_not_assume_all_refs(self):
        ids, _ = polish.complete_reference_ids(SimpleNamespace(), {'core': FastaRecord('core', '', 'M' * 300)})
        self.assertEqual(ids, set())

    def test_complete_metadata_uses_eligibility(self):
        args = SimpleNamespace(reference_metadata={'a': {'reconstruction_allowed': True}, 'b': {'reconstruction_allowed': False}})
        ids, _ = polish.complete_reference_ids(args, {'a': None, 'b': None})
        self.assertEqual(ids, {'a'})

    def test_timeout_repair_budget_and_unique_reference(self):
        a, b = model('a'), model('b', 3000)
        refs = {'core': FastaRecord('core', '', 'M' * 399)}
        args = SimpleNamespace(polish_repair_budget=1, polish_references=2)
        selected, deferred, unavailable = polish.repair_plan([[a], [b]], {'a': {'status': 'timeout'}, 'b': {'status': 'error'}}, refs, {'core'}, args)
        self.assertEqual((len(selected), len(deferred), len(unavailable)), (1, 1, 0))
        self.assertEqual(selected[0]['references'], ['core'])

    def test_zero_budget_records_deferred(self):
        m = model('a', full=False)
        refs = {'core': FastaRecord('core', '', 'M' * 399)}
        selected, deferred, _ = polish.repair_plan([[m]], {}, refs, {'core'}, SimpleNamespace(polish_repair_budget=0))
        self.assertFalse(selected)
        self.assertEqual(len(deferred), 1)

    def test_empty_export_has_all_three_sets_and_header(self):
        with TemporaryDirectory() as directory, patch('orsearcher.engine.validate', return_value={'checked_models': 0}):
            root = Path(directory)
            summary = polish.export_set(root, [], [], SimpleNamespace(sample='insect', genome='unused'))
            self.assertEqual(summary['category_counts'], dict(complete=0, partial=0, pseudogene=0))
            for kind in ('complete', 'partial', 'pseudogene'):
                for suffix in ('proteins.faa', 'cds.fna', 'gff3'):
                    self.assertTrue((root / f'OR_{kind}.{suffix}').exists())
            self.assertIn('candidate_id', (root / 'OR_catalogue.tsv').read_text())

    def test_overlapping_pseudogene_saved_not_counted_independent(self):
        a, p = model('complete'), model('damaged')
        p.update(stops=[1], protein='M' + 'A' * 100 + '*' + 'A' * 297)
        p['_assessment'] = polish.assessment(p, {'core'})
        with TemporaryDirectory() as directory, patch('orsearcher.engine.validate', return_value={'checked_models': 2}):
            summary = polish.export_set(Path(directory), [[a]], [p], SimpleNamespace(sample='insect', genome='unused'))
            self.assertEqual(summary['category_counts']['pseudogene'], 1)
            self.assertEqual(summary['independent_loci'], 1)
            rows = polish.table(Path(directory) / 'OR_catalogue.tsv')
            damage = next(r for r in rows if r['category'] == 'pseudogene')
            self.assertEqual(damage['independent_locus'], 'False')
            self.assertTrue(damage['overlapping_intact_model_ids'])

    def test_load_legacy_pool_without_audit_is_usable(self):
        with TemporaryDirectory() as directory:
            models, pseudo = polish.load_pool(Path(directory), {'g': [model()]}, {}, {'core'})
            self.assertEqual(len(models), 1)
            self.assertFalse(pseudo)

    def test_overlapping_damage_paths_are_one_independent_pseudogene(self):
        p, q, r = model('damaged1'), model('damaged2', 103), model('independent_damage', 5000)
        for m in (p, q, r):
            m['stops'] = [1]
            m['_assessment'] = polish.assessment(m, {'core'})
        with TemporaryDirectory() as directory, patch('orsearcher.engine.validate', return_value={'checked_models': 3}):
            root = Path(directory)
            summary = polish.export_set(root, [], [p, q, r], SimpleNamespace(sample='insect', genome='unused'))
            self.assertEqual(summary['category_counts']['pseudogene'], 3)
            self.assertEqual(summary['independent_loci'], 2)
            rows = polish.table(root / 'OR_catalogue.tsv')
            held = [r for r in rows if r['independent_locus'] == 'False']
            self.assertEqual(len(held), 1)
            self.assertTrue(held[0]['overlapping_pseudogene_model_ids'])

    def test_missing_previous_model_store_is_explicit_error(self):
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, 'all_models.json'):
                polish.load_pool(Path(directory), None, None, {'core'})

    def test_cached_background_is_not_used_when_current_check_disabled(self):
        class FakeEvidence:
            def __init__(self, args):
                self.args = args
            def apply(self, models, label):
                pass
        with TemporaryDirectory() as directory:
            root = Path(directory)
            m = model('a')
            m['background_hits'] = [dict(ref='old_negative', evalue=1e-150, bits=800, identity=95,
                                         qs=1, qe=399, qlen=399, ss=1, se=399, slen=399)]
            args = SimpleNamespace(sample='insect', genome='unused', threads=1, workers=1,
                                   complete_reference_ids={'core'}, polish_repair_budget=0)
            with patch('orsearcher.engine.validate', return_value={'checked_models': 1}):
                summary = polish.polish(args, root, {'core': FastaRecord('core', '', m['protein'])}, FakeEvidence(args), groups={'a': [m]}, decisions={})
            self.assertEqual(summary['category_counts']['complete'], 1)
            self.assertIn('other_protein_background_not_checked', summary['warnings'])
            self.assertTrue((root / 'polish_audit/models.json').exists())
            self.assertIn('background_hits', m)  # Caller/old model was not mutated.

    def test_failed_new_identity_cannot_reuse_cached_positive(self):
        class FakeEvidence:
            def __init__(self, args):
                self.args = args
            def apply(self, models, label):
                for m in models:
                    m.update(homology_hits=[], or_profile_domains=[], preliminary_identity={'status': 'OR_supported'})
        with TemporaryDirectory() as directory:
            args = SimpleNamespace(sample='insect', genome='unused', threads=1, workers=1,
                                   complete_reference_ids={'core'}, polish_repair_budget=0)
            with patch('orsearcher.engine.validate', return_value={'checked_models': 0}):
                summary = polish.polish(args, Path(directory), {'core': FastaRecord('core', '', 'M' * 399)}, FakeEvidence(args), groups={'a': [model()]})
            self.assertEqual(summary['category_counts']['complete'], 0)
            self.assertEqual(summary['alternative_models'], 1)

    def test_repair_window_clipped_for_plus_and_minus(self):
        for strand in ('+', '-'):
            with self.subTest(strand=strand), TemporaryDirectory() as directory:
                root = Path(directory)
                genome = root / 'genome.fna'
                genome.write_text('>chr1\n' + 'A' * 5000 + '\n')
                m = model()
                m['strand'] = strand
                args = SimpleNamespace(genome=genome, max_intron=1000, max_window_bp=1400,
                                       polish_flank_bp=80000, polish_timeout=180)
                job = dict(model=m, candidate_id='case1', references=['core'], original_group_ids=['a'])
                with patch('orsearcher.support.review_or_catalogue.model_candidates', return_value=[]) as rebuild:
                    values, status = polish.repair_case(job, root, args, {'core': FastaRecord('core', '', m['protein'])}, {'chr1': 5000})
                parameters = rebuild.call_args.args[6]
                self.assertEqual(parameters.flank_bp, 100)
                self.assertEqual(status['region_flank_bp'], 100)
                self.assertEqual(values, [])

    def test_exonerate_timeout_retains_original_and_never_creates_pseudogene(self):
        from orsearcher.fasta import write_fasta
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / 'genome.fna'
            genome.write_text('>chr1\n' + 'A' * 5000 + '\n')
            m = model()
            args = SimpleNamespace(genome=genome, max_intron=1000, max_window_bp=1400,
                                   polish_flank_bp=80000, polish_timeout=180)
            job = dict(model=m, candidate_id='case1', references=['core'], original_group_ids=['a'])
            def prepare_region(row, case, *unused, **kwargs):
                write_fasta(case / 'region.fna', [FastaRecord('region', '', 'A' * 1400)])
                polish.save(case / 'region.json', dict(offset=0))
                return []
            with patch('orsearcher.support.review_or_catalogue.model_candidates', side_effect=prepare_region), patch(
                    'orsearcher.support.review_or_catalogue.command', side_effect=subprocess.TimeoutExpired('exonerate', 180)):
                extras, status = polish.repair_case(job, root, args, {'core': FastaRecord('core', '', m['protein'])}, {'chr1': 5000})
            self.assertFalse(extras)
            self.assertEqual(status['issues'][0]['status'], 'timeout')
            clusters, _, _ = polish.consolidate([m, *extras])
            self.assertEqual(len(clusters), 1)
            self.assertTrue(clusters[0][0]['_assessment']['usable_complete'])

    def test_successful_repair_cache_reused_but_timeout_is_retried(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / 'genome.fna'
            genome.write_text('>chr1\n' + 'A' * 5000 + '\n')
            m = model()
            args = SimpleNamespace(genome=genome, max_intron=1000, max_window_bp=1400,
                                   polish_flank_bp=80000, polish_timeout=180)
            job = dict(model=m, candidate_id='case1', references=['core'], original_group_ids=['a'])
            refs = {'core': FastaRecord('core', '', m['protein'])}
            with patch('orsearcher.support.review_or_catalogue.model_candidates', return_value=[]) as rebuild:
                polish.repair_case(job, root, args, refs, {'chr1': 5000})
                _, cached = polish.repair_case(job, root, args, refs, {'chr1': 5000})
                self.assertEqual(rebuild.call_count, 1)
                self.assertTrue(cached['cache_reused'])
                status_path = root / 'repair/cases/case1/status.json'
                state = json.loads(status_path.read_text())
                state['status'] = 'completed_with_warnings'
                polish.save(status_path, state)
                polish.repair_case(job, root, args, refs, {'chr1': 5000})
                self.assertEqual(rebuild.call_count, 2)


if __name__ == '__main__':
    unittest.main()
