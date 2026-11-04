import copy,json,tempfile,unittest,importlib.util,random
from pathlib import Path
from types import SimpleNamespace
from test_identity import protein,hit
from test_v21 import model
from orsearcher.identity import judge
from orsearcher.refine import eligible,plan_groups
from orsearcher.fasta import FastaRecord,write_fasta
from orsearcher.model_store import save_groups,load_groups,pack_groups
from orsearcher.support.audit_or_loci import selection_rank,decide_locus
from orsearcher.support.finalize_catalogue_review import fragment_alias,export_unique


class Candidate08Tests(unittest.TestCase):
    def test_identity_panel_accepts_manifest_string_paths_and_curated_OR_only(self):
        script=Path(__file__).resolve().parents[1]/'scripts/probe_identity_holdout.py'
        spec=importlib.util.spec_from_file_location('identity_probe_test',script)
        probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            for kind,seq in [('OR','M'*150),('GR','C'*150),('other','D'*150)]:
                write_fasta(root/(kind+'.faa'),[FastaRecord(kind,'',seq)])
            result=probe.panel(dict(or_reference=str(root/'OR.faa'),gr_reference=str(root/'GR.faa')),
                {'OR':dict(curated=True)},root/'other.faa',1,random.Random(1))
            self.assertEqual(set(result),{'OR','GR','other'})
            self.assertEqual(result['OR'][0].id,'OR')

    def test_far_OR_supported_by_significance_profile_not_identity_alone(self):
        m=protein();m['homology_hits']=[hit(identity=22.2,bits=104)]
        self.assertEqual(judge(m)['status'],'OR_supported')
        self.assertTrue(judge(m)['signals']['distant_OR_path'])

    def test_far_homolog_without_profile_stays_unknown(self):
        m=protein();m['or_profile_domains']=[];m['homology_hits']=[hit(identity=22.2,bits=104)]
        self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_weak_single_homolog_and_profile_stays_unknown(self):
        m=protein();m['homology_hits']=[hit(bits=40,identity=22,coverage=.4)]
        self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_weak_curated_consensus_and_strong_profile_can_support_OR(self):
        m=protein();m['homology_hits']=[hit(bits=42,identity=24,coverage=.5),hit('or2',bits=40,identity=23,coverage=.5)]
        self.assertEqual(judge(m)['status'],'OR_supported')

    def test_distant_GR_competitor_is_not_ignored(self):
        m=protein();m['homology_hits']=[hit(identity=22,bits=100)]
        m['gr_homology_hits']=[hit('gr',identity=22,bits=160)]
        self.assertNotEqual(judge(m)['status'],'OR_supported')

    def test_qualified_exact_OR_not_vetoed_by_tree_alone(self):
        m=protein();m['qualified_reference_exact_matches']=['or1']
        self.assertEqual(judge(m,{'identity_tree_vote':'GR'})['status'],'OR_supported')

    def test_actual_GR_competition_remains_unresolved(self):
        m=protein();m['gr_homology_hits']=[hit('gr',bits=195)]
        self.assertEqual(judge(m,{'identity_tree_vote':'GR'})['status'],'identity_unresolved')

    def test_native_qualified_alignment_can_enter_repair_without_OR_identity(self):
        m=protein();m.update(or_evalue=1,homology_hits=[],or_profile_domains=[],
            alignment_reference_qualified=True,ref_cov=.7,identity=.2)
        self.assertTrue(eligible(m));self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_scheduling_spreads_budget_before_revisiting_same_region(self):
        a=protein();a.update(contig='chr',ref='r1',exons=[[1,903]])
        b=copy.deepcopy(a);b['exons']=[[4,906]]
        c=copy.deepcopy(a);c['exons']=[[3000,3902]]
        args=SimpleNamespace(max_window_bp=500000,structure_references=3,structure_max_loci=2,
            reference_metadata={'r1':dict(reconstruction_allowed=True)})
        selected,deferred,_=plan_groups([a,b,c],args,{'r1':FastaRecord('r1','','M'*300)})
        self.assertEqual(selected,[a,c]);self.assertEqual(deferred,[b])

    def test_endpoints_precede_small_protein_score_gain(self):
        baseline=model();baseline['terminal']={};alternative=copy.deepcopy(baseline)
        alternative.update(id='different',source='miniprot_terminal',or_bits=205,
            terminal=dict(n_difference=3,c_difference=0))
        self.assertGreater(selection_rank(baseline,[baseline,alternative]),selection_rank(alternative,[baseline,alternative]))
        self.assertIs(decide_locus({},[baseline,alternative],None)['model'],baseline)

    def test_shared_store_roundtrip_preserves_distinct_genomic_placements(self):
        a=protein();a['contig']='chr';b=copy.deepcopy(a);b.update(id='other',cid='other',exons=[[3001,3903]])
        groups={'c':[a],'other':[b]}
        packed=pack_groups(groups);self.assertEqual(len(packed['payloads']),1)
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'models.json';save_groups(path,groups);restored=load_groups(path)
            self.assertEqual(restored,groups)
            restored['c'][0]['homology_hits'][0]['bits']=0
            self.assertEqual(restored['other'][0]['homology_hits'][0]['bits'],200)

    def test_old_model_group_files_still_load(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'old.json';p.write_text(json.dumps({'c':[model()]}))
            self.assertEqual(load_groups(p),{'c':[model()]})

    def test_near_contained_fragment_can_alias_complete_without_length_maximization(self):
        a=model();b=model('miniprot');b.update(exons=[[34,930]],has_start=False,has_stop=False)
        self.assertTrue(fragment_alias(a,b))
        with tempfile.TemporaryDirectory() as d:
            rows=[dict(candidate_id=k,sequence_id='chr') for k in ['a','b']]
            decisions={'a':dict(category='complete',confidence='reviewed_complete',model=a),
                'b':dict(category='partial',confidence='unresolved',model=b)}
            s=export_unique(rows,decisions,Path(d),allow_containment=True)
            self.assertEqual(s['representative_model_records'],1)
            self.assertIn('same_frame_near_contained_complete_fragment',(Path(d)/'nonredundant/OR_candidate_aliases.tsv').read_text())

    def test_damaged_fragment_not_hidden_by_complete_alias(self):
        a=model();b=model('miniprot');b.update(exons=[[34,930]],frameshifts=1)
        self.assertFalse(fragment_alias(a,b))

    def test_incompatible_reading_frame_not_hidden(self):
        a=model();b=model('miniprot');b['exons']=[[35,931]]
        self.assertFalse(fragment_alias(a,b))
