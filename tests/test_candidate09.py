import copy,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from test_identity import protein,hit
from test_v21 import model
import test_v21
from orsearcher.identity import judge,structural_family,full_curated_consensus
from orsearcher.refine import eligible,plan_groups
from orsearcher.support import review_or_catalogue as review,audit_or_loci as audit,crosscheck_or_catalogue as cross
from orsearcher.structure_rescue import reference_choices,rescue_priority,reconstruct,legacy_priority
from orsearcher.fasta import FastaRecord
from orsearcher.cli import parser


def profile_negative():
    m=protein();m['or_profile_domains']=[]
    m['homology_hits']=[hit('whole',bits=370,identity=44,coverage=.98),hit('secondary',bits=57,identity=26,coverage=.4)]
    return m


class Candidate09Tests(unittest.TestCase):
    def test_whole_curated_plus_distinct_homolog_stable_to_tree_votes(self):
        for vote in ['OR','GR','unresolved']:
            r=judge(profile_negative(),dict(identity_tree_vote=vote))
            self.assertEqual(r['status'],'OR_supported')
            self.assertTrue(r['signals']['profile_independent_full_consensus'])

    def test_one_strong_whole_hit_remains_insufficient(self):
        m=profile_negative();m['homology_hits']=m['homology_hits'][:1]
        self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_secondary_alias_is_not_independent_consensus(self):
        m=profile_negative();m['homology_hits'][1]['reference_sequence_key']='whole'
        self.assertFalse(full_curated_consensus(m['homology_hits']))

    def test_whole_hit_must_cover_both_query_and_subject(self):
        for field in ['qe','se']:
            m=profile_negative();m['homology_hits'][0][field]=150
            self.assertFalse(judge(m)['signals']['profile_independent_full_consensus'])

    def test_low_identity_full_hit_not_profile_independent(self):
        m=profile_negative();m['homology_hits'][0]['identity']=25
        self.assertFalse(full_curated_consensus(m['homology_hits']))

    def test_unreviewed_whole_reference_not_curated_consensus(self):
        m=profile_negative();m['homology_hits'][0]['curated']=False
        self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_GR_profile_blocks_new_profile_negative_path(self):
        m=profile_negative();m['gr_profile_domains']=[[dict(domain_evalue=1e-20,coverage=.8)]]
        self.assertFalse(judge(m)['signals']['profile_independent_full_consensus'])
        self.assertNotEqual(judge(m)['status'],'OR_supported')

    def test_close_GR_competition_blocks_new_path(self):
        m=profile_negative();m['gr_homology_hits']=[hit('gr',bits=355)]
        self.assertFalse(judge(m)['signals']['profile_independent_full_consensus'])
        self.assertNotEqual(judge(m)['status'],'OR_supported')

    def test_background_competition_not_silently_overridden(self):
        m=profile_negative();m['background_hits']=[hit('negative',bits=600,identity=95)]
        self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_profile_negative_identity_does_not_ignore_coding_damage(self):
        m=profile_negative();m.update(or_ga=False,or_evalue=1,full_ref_support=['whole'])
        structural_family(m)
        self.assertTrue(review.complete_model(m))
        m['frameshifts']=1;self.assertFalse(review.complete_model(m))

    def test_structural_eligibility_not_final_identity_confirmation(self):
        m=protein();m['homology_hits']=[hit(bits=70,identity=20)];m['or_profile_domains']=[]
        structural_family(m)
        self.assertTrue(m['structure_or_eligible']);self.assertEqual(m['family'],'OR_GR_uncertain')

    def test_weak_curated_homolog_can_enter_repair_not_identity(self):
        m=protein();m.update(or_evalue=1,or_profile_domains=[])
        m['homology_hits']=[hit(bits=32,identity=22,coverage=.7)]
        self.assertTrue(eligible(m));self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_partial_rank_compares_same_reference_not_native_identity(self):
        short=model();long=copy.deepcopy(short)
        short.update(has_start=False,has_stop=False,length=150,id='short',identity=.95,ref_cov=.9,full_ref_support=[],homology_hits=[hit('shared',bits=100,coverage=.5)])
        long.update(has_start=False,has_stop=False,length=210,id='long',identity=.25,ref_cov=.4,full_ref_support=[],homology_hits=[hit('shared',bits=145,coverage=.7)])
        self.assertGreater(audit.partial_selection_rank(long,[short,long]),audit.partial_selection_rank(short,[short,long]))
        self.assertIs(audit.decide_locus({},[dict(short,source='original'),long],None)['model'],long)

    def test_length_without_extra_homolog_support_has_no_rank_effect(self):
        a=model();b=copy.deepcopy(a);b['length']=600
        self.assertEqual(audit.partial_selection_rank(a,[a,b]),audit.partial_selection_rank(b,[a,b]))

    def test_noncanonical_longer_model_not_forced_into_representative(self):
        a=model();a.update(has_start=False,full_ref_support=[])
        b=copy.deepcopy(a);b.update(source='exonerate',id='long',canonical=False,length=500,or_bits=900)
        self.assertIs(audit.decide_locus({},[a,b],None)['model'],a)

    def test_damaged_high_score_path_not_selected_over_intact_partial(self):
        a=model();a.update(has_start=False,full_ref_support=[])
        b=copy.deepcopy(a);b.update(source='exonerate',id='damaged',frameshifts=1,or_bits=900)
        self.assertIs(audit.decide_locus({},[a,b],None)['model'],a)

    def test_strong_GR_not_structurally_eligible_as_OR(self):
        m=protein();m['gr_homology_hits']=[hit('gr',bits=600)]
        m['gr_profile_domains']=[[dict(domain_evalue=1e-20,coverage=.8)]]
        structural_family(m);self.assertFalse(m['structure_or_eligible'])

    def test_rescue_references_distinct_audited_and_unused(self):
        m=protein();m['homology_hits']=[hit('used'),hit('alias'),hit('new'),hit('bad')]
        refs={r:FastaRecord(r,'',s*300) for r,s in [('used','A'),('alias','A'),('new','C'),('bad','D')]}
        meta={r:dict(reconstruction_allowed=r!='bad') for r in refs}
        self.assertEqual(list(reference_choices([m],refs,meta,{'A'*300},3)),['new'])

    def test_rescue_reference_budget_respects_sequence_duplicates(self):
        m=protein();m['homology_hits']=[hit('a'),hit('alias'),hit('b')]
        refs={r:FastaRecord(r,'',s*300) for r,s in [('a','A'),('alias','A'),('b','C')]}
        meta={r:dict(reconstruction_allowed=True) for r in refs}
        self.assertEqual(len(reference_choices([m],refs,meta,set(),1)),1)
        self.assertEqual(len(reference_choices([m],refs,meta,set(),3)),2)

    def test_complete_model_does_not_consume_extra_reconstruction_budget(self):
        self.assertIsNone(rescue_priority([protein()]))

    def test_legacy_budget_prioritizes_incomplete_over_high_score_complete(self):
        full=protein();full['or_bits']=1000
        partial=protein();partial.update(has_stop=False,or_bits=100)
        self.assertGreater(legacy_priority([partial]),legacy_priority([full]))

    def test_close_disruption_still_gets_complementary_priority(self):
        disrupted=protein();disrupted.update(frameshifts=1,ref_cov=.99,identity=.99)
        self.assertGreater(legacy_priority([disrupted]),legacy_priority([protein()]))

    def test_credible_incomplete_model_can_request_rescue(self):
        m=protein();m.update(has_start=False,contig='chr');self.assertIsNotNone(rescue_priority([m]))

    def test_disabled_rescue_does_not_run_tools(self):
        with tempfile.TemporaryDirectory() as d,patch('orsearcher.structure_rescue.cross.process') as run:
            result=reconstruct(SimpleNamespace(structure_rescue_budget=0),[],{},Path(d),None,{})
            self.assertEqual(result['status'],'disabled');run.assert_not_called()

    def test_primary_segfault_retry_uses_remaining_timeout(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=test_v21.V21Tests().fixture(root);calls=[]
            (root/'run/cases/c/miniprot.gff3').write_text('##gff-version 3\n')
            def command(cmd,output,timeout):
                calls.append(timeout)
                if len(calls)==1:raise RuntimeError('exit=-11')
                Path(output).write_text('')
            args=SimpleNamespace(references_per_candidate=1,pseudogene_references=1,max_window_bp=2000,
                exonerate='exonerate',max_intron=1000,timeout=20,legacy_rescue='off')
            with patch('orsearcher.support.review_or_catalogue.command',command):
                _,status=cross.process(row,[model()],root/'run',root/'out',dict(selected_model_id=model()['id'],model_category='partial'),args)
            self.assertEqual(status['status'],'completed');self.assertLessEqual(calls[1],calls[0])

    def test_rescue_stays_same_aligner(self):
        self.assertEqual(audit.tool(model('exonerate_rescue')),'exonerate')

    def test_new_aligner_models_have_fusion_flags_before_rescue(self):
        from orsearcher.engine import initialize_fusion_evidence
        m=protein();m.update(has_start=False,contig='chr');del m['fusion_pairs']
        with patch('orsearcher.engine.audit.anchor_groups',return_value=[]),patch('orsearcher.engine.audit.fusion_links',return_value=[]):
            initialize_fusion_evidence([m])
        self.assertEqual(m['fusion_pairs'],[]);self.assertIsNotNone(rescue_priority([m]))

    def test_cli_exposes_bounded_rescue_controls(self):
        args=parser().parse_args(['refine','--genome','g','--outdir','o','--or-reference','r','--gr-reference','s',
            '--or-hmm','h','--gr-hmm','i','--candidate-models','j'])
        self.assertEqual(args.structure_rescue_budget,64);self.assertEqual(args.structure_rescue_timeout,20)

    def test_rescue_budget_deferral_and_sources_are_explicit(self):
        from orsearcher.fasta import write_fasta
        from orsearcher.engine import row_for
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);byid={};rows=[]
            for cid in ['a','b','c']:
                m=protein();m.update(cid=cid,id=cid+'.original',contig='chr',has_start=False)
                byid[cid]=[m];rows.append(row_for(m));case=root/'cases'/cid;case.mkdir(parents=True)
                write_fasta(case/'region.fna',[FastaRecord('region','','A'*1000)])
                write_fasta(case/'references.faa',[FastaRecord('unused','','M'*300)])
                (case/'region.json').write_text('{}');(case/'miniprot.gff3').write_text('##gff-version 3\n')
            args=SimpleNamespace(structure_rescue_budget=1,structure_rescue_references=2,structure_rescue_timeout=20,
                exonerate_timeout=30,workers=1,max_intron=1000,max_window_bp=500000,
                reference_metadata={'or1':dict(reconstruction_allowed=True)})
            def process(row,models,*rest):
                return [dict(models[0],source='exonerate',id=row['candidate_id']+'.extra')],dict(status='completed',candidate_id=row['candidate_id'])
            evidence=SimpleNamespace(apply=lambda *rest:None)
            with patch('orsearcher.structure_rescue.cross.process',side_effect=process) as run:
                result=reconstruct(args,rows,byid,root,evidence,{'or1':FastaRecord('or1','','C'*300)})
            self.assertEqual(run.call_count,1);self.assertEqual(result['selected'],1);self.assertEqual(result['deferred'],2)
            extra=[m for ms in byid.values() for m in ms if m['source']=='exonerate_rescue']
            self.assertEqual(len(extra),1);self.assertEqual(audit.tool(extra[0]),'exonerate')
            self.assertTrue((root/'rescue/plan.json').exists())
