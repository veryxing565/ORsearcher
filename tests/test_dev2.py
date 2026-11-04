import copy,json,tempfile,unittest,subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from test_v21 import model
import test_v21
from orsearcher.evidence import alignment_provenance
from orsearcher.support import audit_or_loci as audit,crosscheck_or_catalogue as cross
from orsearcher.support.finalize_catalogue_review import fragment_alias,export_unique

class Dev2Tests(unittest.TestCase):
    def test_self_and_annotation_alignment_are_not_external(self):
        for ref in ['adaptive_1','annotation_1']:
            m=model(ref=ref);alignment_provenance(m,{})
            self.assertFalse(audit.external_alignment(m));self.assertEqual(audit.external_identity(m),0)
        m=model();alignment_provenance(m,{'r1':dict(reconstruction_allowed=True,sequence_key='seq1')})
        self.assertTrue(audit.external_alignment(m));self.assertEqual(audit.reference_key(m),'seq1')
    def test_self_match_cannot_upgrade_complete_confidence(self):
        a=model(ref='adaptive_1');b=model('exonerate','r2');a['identity']=1.
        d=audit.decide_locus({},[a,b],None)
        self.assertEqual(d['confidence'],'provisional_complete');self.assertIs(d['model'],b)
    def test_votes_deduplicate_reference_sequence_and_terminal_variants(self):
        a=model();bs=[model('exonerate','r'+str(i)) for i in range(2,8)]
        for b in bs:b['reference_sequence_key']='same_external_sequence'
        self.assertEqual(len(audit.compatible_support(a,bs)),1)
        self.assertEqual(audit.selection_rank(a,[a]+bs)[2],1)
    def test_external_exact_protein_beats_longer_self_seed(self):
        a=model('exonerate');a['reference_exact_matches']=['r1'];a['or_bits']=798
        b=model('miniprot','adaptive_49');b.update(exons=[[1,957]],length=318,identity=1.,or_bits=787)
        d=audit.decide_locus({},[dict(a,source='original'),a,b],None)
        self.assertEqual(d['model']['length'],300);self.assertEqual(d['confidence'],'provisional_complete')
    def test_full_external_alignment_beats_long_terminal_inference(self):
        a=model();a.update(ref_cov=.89,identity=.27,or_bits=129)
        b=model('exonerate','r2');b.update(exons=[[31,903]],length=290,ref_cov=.56,identity=.30,or_bits=129)
        d=audit.decide_locus({},[a,b],None)
        self.assertIs(d['model'],a)
    def test_length_is_not_a_selection_tie_breaker(self):
        a=model();b=copy.deepcopy(a);b['length']=340
        self.assertEqual(audit.selection_rank(a,[]),audit.selection_rank(b,[]))
    def test_more_protein_hits_do_not_overrule_better_external_model(self):
        a=model();a.update(ref_cov=.987,identity=.2642,or_bits=159,full_ref_support=['r1'])
        b=model();b.update(exons=[[37,903]],length=288,ref_cov=.862,identity=.2663,or_bits=157,full_ref_support=['r1','r2'])
        d=audit.decide_locus({},[a,b],None)
        self.assertIs(d['model'],a)
    def test_comparable_internal_splice_paths_are_unresolved(self):
        a=model();a.update(exons=[[1,300],[401,1003]],identity=.35)
        b=model('exonerate','r2');b.update(exons=[[1,303],[404,1003]],identity=.35)
        d=audit.decide_locus({},[a,b],None)
        self.assertEqual(d['category'],'partial');self.assertIn('internal_splice',d['reason'])
    def test_exact_external_splice_path_can_remain_complete(self):
        a=model();a.update(exons=[[1,300],[401,1003]],identity=.35,reference_exact_matches=['r1'])
        b=model('exonerate','r2');b.update(exons=[[1,303],[404,1003]],identity=.35)
        d=audit.decide_locus({},[a,b],None)
        self.assertEqual(d['category'],'complete');self.assertIs(d['model'],a)
    def test_partial_alias_requires_same_frame_and_intact_splicing(self):
        a=model();b=copy.deepcopy(a);b.update(exons=[[304,903]],cds=a['cds'][303:],protein=a['protein'][101:],length=199,has_start=False)
        self.assertTrue(fragment_alias(a,b))
        b['exons']=[[305,903]];self.assertFalse(fragment_alias(a,b))
        b['exons']=[[304,903]];b['frameshifts']=1;self.assertFalse(fragment_alias(a,b))
    def test_fragment_alias_does_not_bridge_two_tandem_copies(self):
        a=model();b=copy.deepcopy(a);b.update(exons=[[1,1806]],length=601)
        c=copy.deepcopy(a);c['exons']=[[904,1806]]
        self.assertFalse(fragment_alias(a,b));self.assertFalse(fragment_alias(c,b));self.assertFalse(fragment_alias(a,c))
    def test_fragment_alias_rejects_intron_retention(self):
        a=model();a['exons']=[[1,300],[355,903]]
        b=copy.deepcopy(a);b['exons']=[[1,903]]
        self.assertFalse(fragment_alias(a,b))
    def test_partial_alias_is_documented_and_raw_models_retained(self):
        a=model();b=copy.deepcopy(a);b.update(exons=[[364,903]],cds=a['cds'][363:],protein=a['protein'][121:],length=179,has_start=False,or_bits=120)
        rows=[dict(candidate_id=k,sequence_id='chr') for k in ['a','b']]
        ds={k:dict(category='partial',confidence='unresolved',model=m) for k,m in [('a',a),('b',b)]}
        with tempfile.TemporaryDirectory() as d:
            stats=export_unique(rows,ds,Path(d),allow_containment=True)
            self.assertEqual(stats['representative_loci'],1);self.assertIn('same_frame_near_contained_fragment',(Path(d)/'nonredundant/OR_candidate_aliases.tsv').read_text())
            self.assertEqual(len(ds),2)
    def test_complementary_seed_is_still_the_same_tool(self):
        self.assertEqual(audit.tool(model('exonerate_legacy','r2')),'exonerate')
    def test_duplicate_query_proteins_are_one_reconstruction_reference(self):
        from orsearcher.fasta import FastaRecord
        refs={k:FastaRecord(k,'',seq) for k,seq in [('r1','MAAAA'),('r2','MAAAA'),('r3','MCCCC')]}
        a=model();_,names=cross.target_references([a],refs,dict(selected_model_id=a['id']),2)
        self.assertEqual(names,['r1','r3'])
    def test_strong_full_homolog_not_omitted_by_short_identity(self):
        from orsearcher.fasta import FastaRecord
        refs={k:FastaRecord(k,'',seq) for k,seq in [('local','MAAAA'),('whole','MCCCC'),('other','MGGGG')]}
        a=model(ref='local');a.update(ref_cov=.5,identity=.95,full_ref_support=['whole'])
        _,names=cross.target_references([a],refs,dict(selected_model_id=a['id'],model_category='complete'),2)
        self.assertEqual(names,['whole','local'])
    def test_legacy_segfault_retries_within_remaining_limit(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=test_v21.V21Tests().fixture(root);(root/'run/cases/c/miniprot.gff3').write_text('##gff-version 3\n');calls=[]
            def command(cmd,output,timeout):
                calls.append((cmd,timeout))
                if cmd[cmd.index('--proteinwordlen')+1]=='6' and cmd[cmd.index('--refine')+1]=='region':raise RuntimeError('exit=-11')
                Path(output).write_text('')
            args=SimpleNamespace(references_per_candidate=1,pseudogene_references=1,max_window_bp=2000,exonerate='exonerate',max_intron=1000,timeout=60,legacy_rescue='auto',exonerate_word_length=3)
            with patch('orsearcher.support.review_or_catalogue.command',command):
                _,status=cross.process(row,[model()],root/'run',root/'out',dict(selected_model_id=model()['id'],model_category='complete'),args)
            self.assertEqual(len(calls),3);self.assertLessEqual(calls[-1][1],30)
            self.assertEqual(calls[-1][0][calls[-1][0].index('--refine')+1],'none')
            self.assertTrue(status['legacy_fallback_without_refinement']);self.assertEqual(status['status'],'completed')
    def test_legacy_rescue_is_bounded_and_both_outputs_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=test_v21.V21Tests().fixture(root);(root/'run/cases/c/miniprot.gff3').write_text('##gff-version 3\n')
            calls=[]
            def command(cmd,output,timeout):calls.append((cmd,timeout));Path(output).write_text('')
            args=SimpleNamespace(references_per_candidate=1,pseudogene_references=1,max_window_bp=2000,exonerate='exonerate',max_intron=1000,timeout=60,legacy_rescue='auto',exonerate_word_length=3)
            with patch('orsearcher.support.review_or_catalogue.command',command):
                _,status=cross.process(row,[model()],root/'run',root/'out',dict(selected_model_id=model()['id'],model_category='complete'),args)
            self.assertEqual(len(calls),2);self.assertEqual(calls[1][1],30)
            self.assertEqual(calls[1][0][calls[1][0].index('--proteinwordlen')+1],'6')
            self.assertEqual(status['legacy_status'],'completed')
            self.assertTrue((root/'out/cases/c/exonerate.out').exists());self.assertTrue((root/'out/cases/c/exonerate_legacy.out').exists())
    def test_legacy_budget_deferral_is_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=test_v21.V21Tests().fixture(root);(root/'run/cases/c/miniprot.gff3').write_text('##gff-version 3\n')
            def command(cmd,output,timeout):Path(output).write_text('')
            args=SimpleNamespace(references_per_candidate=1,pseudogene_references=1,max_window_bp=2000,exonerate='exonerate',max_intron=1000,timeout=60,legacy_rescue='auto',exonerate_word_length=3,allow_legacy=False)
            with patch('orsearcher.support.review_or_catalogue.command',command):
                _,status=cross.process(row,[model()],root/'run',root/'out',dict(selected_model_id=model()['id'],model_category='complete'),args)
            self.assertEqual(status['legacy_status'],'deferred_budget');self.assertEqual(status['status'],'completed')
