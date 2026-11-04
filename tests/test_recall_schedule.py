import copy,tempfile,unittest,json
from pathlib import Path
from types import SimpleNamespace
from orsearcher.fasta import FastaRecord
from orsearcher.sharding import Segment
from orsearcher.recall import unique_anchors,local_hints,oversized_hints,substantial,rebuild
from unittest.mock import patch
from orsearcher.recall_schedule import plan,exon_blocks,region,reference_names,origin,spread_order
from orsearcher.translated import reconstruct


def hint(s,e=None,source='metaeuk',score=100,ref='r'):
    return dict(contig='chr',strand='+',start=s,end=e or s+99,references=[ref],
                origin=source,bits=score,eligible=True)


class ScheduleTests(unittest.TestCase):
    def setup_args(self,budget=2):
        return SimpleNamespace(recall_flank_bp=5000,recall_window_bp=60000,recall_max_regions=budget,
            recall_references=4,reference_metadata={'r':{'reconstruction_allowed':True}},max_intron=200000)
    def test_long_anchor_cannot_swallow_contained_short_anchor(self):
        a=hint(1,300000,score=1000);b=hint(1000,1300,source='translated_hsp',score=60)
        self.assertEqual(len(unique_anchors([a,b])),2)
    def test_reciprocal_overlap_requires_both_spans(self):
        a=hint(1,1000);b=hint(1,100)
        self.assertEqual(len(unique_anchors([a,b])),2)
    def test_identical_cross_source_hints_keep_source_provenance(self):
        a=hint(1);b=hint(1,source='translated_hsp',score=30)
        result=unique_anchors([a,b]);self.assertEqual(len(result),1)
        self.assertEqual(set(result[0]['origins']),{'metaeuk','translated_hsp'})
        self.assertEqual(origin(result[0]),'translated_hsp')
    def test_oversized_regions_do_not_consume_budget(self):
        a=self.setup_args();refs={'r':FastaRecord('r','','M'*300)}
        seeds=[hint(1,300000,score=1000),hint(400000),hint(500000)]
        selected,deferred,skipped,summary=plan(a,seeds,refs,{'chr':1000000})
        self.assertEqual(len(selected),2);self.assertFalse(deferred)
        self.assertEqual(summary['skipped_counts'],{'window_too_large':1})
    def test_missing_references_do_not_consume_budget(self):
        a=self.setup_args(1);refs={'r':FastaRecord('r','','M'*300)}
        selected,_,skipped,_=plan(a,[hint(1,ref='absent'),hint(20000)],refs,{'chr':1000000})
        self.assertEqual(len(selected),1);self.assertEqual(skipped[0]['status'],'no_reference')
    def test_source_round_robin_not_mixed_native_score(self):
        a=self.setup_args();refs={'r':FastaRecord('r','','M'*300)}
        seeds=[hint(1,score=9999),hint(20000,score=9998),hint(40000,source='translated_hsp',score=30)]
        selected,_,_,summary=plan(a,seeds,refs,{'chr':1000000})
        self.assertEqual(summary['selected_origins'],{'translated_hsp':1,'metaeuk':1})
    def test_spatial_rounds_do_not_drop_dense_tandem_hints(self):
        seeds=[hint(1,score=100),hint(1000,score=99),hint(20000,score=98)]
        result=spread_order(seeds);self.assertEqual([a['start'] for a in result],[1,20000,1000])
        self.assertEqual(len(result),3)
    def test_blocks_split_at_exon_boundaries_and_preserve_all_exons(self):
        exons=[[1,300],[1000,1299],[200000,200299]]
        blocks=exon_blocks(exons,50000);self.assertEqual(len(blocks),2)
        self.assertEqual([e for b in blocks for e in b],list(map(tuple,exons)))
    def test_individual_huge_exon_not_silently_truncated(self):
        self.assertEqual(exon_blocks([[1,100000]],50000),[[(1,100000)]])
    def test_chain_fallback_to_local_window(self):
        a=hint(200000);a.update(chain_start=1,chain_end=400000)
        self.assertEqual(region(a,1000000,5000,60000),(195000,205099))
    def test_contig_edge_clamping(self):
        self.assertEqual(region(hint(1),1000,5000,60000),(1,1000))
    def test_complementary_uses_checked_local_fallback_not_long_chain(self):
        a=hint(200000);a.update(chain_start=1,chain_end=400000,planned_start=195000,planned_end=205099)
        args=SimpleNamespace(genome=Path('genome.fna'),translated_flank_bp=5000,max_window_bp=60000,reference_metadata={})
        response=SimpleNamespace(stdout='>region\n'+'A'*10100+'\n')
        with tempfile.TemporaryDirectory() as d,patch('orsearcher.translated.subprocess.run',return_value=response) as run:
            _,status=reconstruct(args,a,{},Path(d),{'chr':1000000})
            self.assertEqual(run.call_args.args[0][-1],'chr:195000-205099')
            self.assertEqual(status['status'],'no_reconstruction_reference')
    def test_complementary_uses_clamped_short_contig_window(self):
        a=hint(1);a.update(planned_start=1,planned_end=1000)
        args=SimpleNamespace(genome=Path('genome.fna'),translated_flank_bp=5000,max_window_bp=60000,reference_metadata={})
        response=SimpleNamespace(stdout='>region\n'+'A'*1000+'\n')
        with tempfile.TemporaryDirectory() as d,patch('orsearcher.translated.subprocess.run',return_value=response) as run:
            _,status=reconstruct(args,a,{},Path(d),{'chr':1000})
            self.assertEqual(run.call_args.args[0][-1],'chr:1-1000')
            self.assertEqual(status['status'],'no_reconstruction_reference')
    def test_invalid_planned_window_fails_before_external_command(self):
        a=hint(1);a.update(planned_start=0,planned_end=1000)
        args=SimpleNamespace(genome=Path('genome.fna'),translated_flank_bp=5000,max_window_bp=60000)
        with tempfile.TemporaryDirectory() as d,patch('orsearcher.translated.subprocess.run') as run:
            with self.assertRaises(ValueError):reconstruct(args,a,{},Path(d),{'chr':1000})
            run.assert_not_called()
    def test_distinct_reference_proteins_only(self):
        a=hint(1);a['references']=['r','s','t'];refs={k:FastaRecord(k,'',seq) for k,seq in [('r','MA'),('s','MA'),('t','MC')]}
        meta={k:{'reconstruction_allowed':True} for k in refs}
        self.assertEqual(reference_names(a,refs,meta,4),['r','t'])
    def test_split_parent_profile_not_broadcast_as_local_evidence(self):
        a=self.setup_args();m=dict(id='large',contig='chr',strand='+',exons=[[1,300],[100000,100299]],
            ref='r',native_score=100,or_bits=300,or_ga=True,or_evalue=1e-30)
        hints=local_hints(a,[m],[],[])
        self.assertEqual(len(hints),2);self.assertTrue(all(not h['or_signal'] and h['or_bits']==0 for h in hints))
        self.assertEqual(m['exons'],[[1,300],[100000,100299]])
    def test_oversized_raw_call_maps_back_to_original_contig(self):
        seg=Segment('v','chr','',1000,501000,1000000,'s',True)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'s';p.mkdir();(p/'skipped_models.json').write_text(json.dumps([{'model':dict(virtual_id='v',strand='-',ref='R',native_score=40,exons=[[1,300],[100000,100299]])}]))
            h=oversized_hints(Path(d),{'v':seg},{'R':'r'})[0]
            self.assertEqual(h['exons'],[[1001,1300],[101000,101299]])
            self.assertTrue(h['hint_only']);self.assertNotIn('cds',h)
    def test_plan_does_not_mutate_input_references(self):
        a=self.setup_args();seed=hint(1);before=copy.deepcopy(seed)
        plan(a,[seed],{'r':FastaRecord('r','','M'*300)},{'chr':1000000})
        self.assertEqual(seed,before)
    def test_partial_model_cannot_skip_complementary_reconstruction(self):
        m=dict(ref_cov=.95,length=350,has_start=True,has_stop=False,canonical=True)
        self.assertFalse(substantial(m))
    def test_complete_looking_model_only_skips_scheduling(self):
        m=dict(ref_cov=.95,length=350,has_start=True,has_stop=True,canonical=True,frameshifts=0,stops=[],ambiguous=0,remainder=0)
        self.assertTrue(substantial(m));self.assertNotIn('category',m)
    def test_damaged_model_gets_complementary_attempt(self):
        m=dict(ref_cov=.95,length=350,has_start=True,has_stop=True,canonical=True,frameshifts=1,stops=[],ambiguous=0,remainder=0)
        self.assertFalse(substantial(m))
    def test_complementary_allowance_counts_only_needy_cases(self):
        from test_v21 import model
        a=self.setup_args(3);a.workers=1;a.sample='test';a.recall_complementary_regions=1;a.recall_timeout=20
        a.exonerate_word_length=3
        full=model();full['ref_cov']=.95
        seeds=[hint(1,source='translated_hsp'),hint(20000,source='translated_hsp'),hint(40000,source='translated_hsp')]
        def wise(args,anchor,*rest):return ([full] if anchor['start']==1 else []),dict(status='completed',errors=[])
        def complementary(*rest):return [],dict(status='no_model',errors=[])
        class FakeEvidence:
            def apply(self,*rest):pass
        with tempfile.TemporaryDirectory() as d,patch('orsearcher.recall.local_hints',return_value=seeds),patch('orsearcher.recall.wise_case',side_effect=wise),patch('orsearcher.recall.complementary_reconstruct',side_effect=complementary) as extra:
            _,_,summary=rebuild(a,[],[],[],{'r':FastaRecord('r','','M'*300)},{'chr':1000000},Path(d),FakeEvidence())
            self.assertEqual(summary['complementary']['needed'],2)
            self.assertEqual(summary['complementary']['selected'],1);self.assertEqual(extra.call_count,1)
            self.assertEqual(extra.call_args.args[1]['start'],20000)
