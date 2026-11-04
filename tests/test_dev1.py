import copy,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from test_v21 import model
from orsearcher.fasta import FastaRecord
from orsearcher.translated import anchors,chain_references,compact_models,discover
from orsearcher.sharding import Segment
from orsearcher.discovery import union_models
from orsearcher.engine import initial_case
from orsearcher.support.review_or_catalogue import model_candidates
from orsearcher.support.finalize_catalogue_review import terminal_equivalent,export_unique

def h(s,e,qs,qe,ref='r1',strand='+',bits=70):
    return dict(contig='chr',strand=strand,start=s,end=e,qstart=qs,qend=qe,qlen=400,
                ref=ref,bits=bits,length=qe-qs+1,evalue=1e-8)

class Dev1Tests(unittest.TestCase):
    def test_chain_multiple_exons_both_strands(self):
        for strand in ['+','-']:
            hs=[h(1,300,1 if strand=='+' else 301,100 if strand=='+' else 400,strand=strand),
                h(501,800,301 if strand=='+' else 1,400 if strand=='+' else 100,strand=strand)]
            seed=chain_references(anchors(hs,[],{}),hs,200000,470000)
            self.assertTrue(all(s['reference_chains'][0]['hsp_count']==2 for s in seed))
            self.assertTrue(all(s['reference_chains'][0]['query_coverage']==.5 for s in seed))
    def test_tandem_query_reset_not_a_chain(self):
        hs=[h(1,900,1,300),h(2001,2900,1,300)]
        seed=chain_references(anchors(hs,[],{}),hs,200000,470000)
        self.assertTrue(all(s['reference_chains'][0]['hsp_count']==1 for s in seed))
        self.assertEqual(len(seed),2)
    def test_far_genomic_chain_is_bounded(self):
        hs=[h(1,300,1,100),h(10001,10300,101,200)]
        seed=chain_references(anchors(hs,[],{}),hs,200000,1000)
        self.assertTrue(all(s['reference_chains'][0]['hsp_count']==1 for s in seed))
    def test_chain_upgrades_weak_seed_without_single_hsp_filter(self):
        hs=[h(1,180,1,60,bits=45),h(501,680,61,120,bits=45)]
        seed=anchors(hs,[],{});self.assertFalse(any(s['eligible'] for s in seed))
        chain_references(seed,hs,200000,470000);self.assertTrue(all(s['eligible'] for s in seed))
    def test_original_metadata_precedes_terminal_extension(self):
        # Upstream ATG 51 residues away must NOT be used when alignment starts
        # at residue 1. Older code extended it because ref/alignment were blank.
        seq='ATG'+'GCT'*349+'TAA';exons=[(154,len(seq),0)]
        row=dict(candidate_id='c',sequence_id='chr',start=154,end=len(seq),strand='+',frameshifts=0)
        original=dict(ref='r1',ref_cov=.75,identity=.3,alignment_tool='exonerate',alignment=dict(qstart=1,qend=300,qlen=300))
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);args=SimpleNamespace(flank_bp=1000,samtools='samtools',genome='fake',miniprot='miniprot',threads_per_job=1,max_intron=200000,miniprot_timeout=10)
            with patch('orsearcher.support.review_or_catalogue.command',return_value='>region\n'+seq+'\n'),patch('orsearcher.support.review_or_catalogue.parse_miniprot_gff',return_value=[]):
                values=model_candidates(row,root,exons,'A'*299,{'r1':FastaRecord('r1','','M'+'A'*299)},{'chr':len(seq)},args,original_model=original)
            self.assertEqual(values[0]['ref'],'r1');self.assertEqual(values[0]['alignment_tool'],'exonerate')
            self.assertFalse(any(v['exons'][0][0]==1 for v in values[1:]))
    def test_compact_preserves_distinct_aligner_and_splice_models(self):
        a=model();a.update(contig='chr',alignment_tool='miniprot')
        b=copy.deepcopy(a);b.update(source='exonerate',ref='r2',alignment_tool='exonerate')
        c=copy.deepcopy(a);c['exons']=[[1,300],[401,1003]]
        result=compact_models([a,b,c]);self.assertEqual(len(result),1)
        self.assertEqual(len(result[0]['discovery_alternatives']),2)
        self.assertEqual(compact_models([c,b,a]),result)
    def test_union_does_not_drop_independent_duplicate_fit(self):
        a=model();a.update(contig='chr',search_passes=['first'],search_references=['r1'],alignment_tool='miniprot')
        b=copy.deepcopy(a);b.update(ref='r2',alignment_tool='exonerate',search_passes=['translated'],search_references=['r2'])
        result=union_models([a],[b]);self.assertEqual(len(result),1)
        self.assertEqual(result[0]['discovery_alternatives'][0]['alignment_tool'],'exonerate')
    def test_near_terminal_models_same_frame_both_strands(self):
        for strand in ['+','-']:
            a=model();a['strand']=strand;a['exons']=[[101,1003]]
            b=copy.deepcopy(a);b['exons']=[[95,1003]] if strand=='+' else [[101,1009]]
            self.assertTrue(terminal_equivalent(a,b))
            b['exons']=[[96,1003]] if strand=='+' else [[101,1008]]
            self.assertFalse(terminal_equivalent(a,b))
    def test_terminal_alias_rejects_damage_and_changed_splicing(self):
        a=model();a['exons']=[[1,300],[401,1003]]
        b=copy.deepcopy(a);b['exons']=[[1,303],[404,1003]]
        self.assertFalse(terminal_equivalent(a,b))
        b=copy.deepcopy(a);b['frameshifts']=1;self.assertFalse(terminal_equivalent(a,b))
    def test_near_complete_fragment_export_alias_recorded(self):
        a=model();a['exons']=[[101,1003]]
        b=copy.deepcopy(a);b.update(exons=[[95,1003]],protein='AA'+a['protein'],has_start=False)
        rows=[dict(candidate_id=k,sequence_id='chr') for k in ['a','b']]
        ds={'a':dict(category='complete',confidence='reviewed_complete',model=a),
            'b':dict(category='partial',confidence='unresolved',model=b)}
        with tempfile.TemporaryDirectory() as d:
            stats=export_unique(rows,ds,Path(d),allow_containment=True)
            self.assertEqual(stats['representative_loci'],1)
            self.assertIn('same_frame_terminal_variant',(Path(d)/'nonredundant/OR_candidate_aliases.tsv').read_text())
    def test_discovery_alternatives_reach_structural_review(self):
        a=model();a.update(contig='chr',alignment_tool='miniprot')
        b=copy.deepcopy(a);b.update(ref='r2',alignment_tool='exonerate',source='exonerate')
        a['discovery_alternatives']=[b]
        row=dict(candidate_id='c',sequence_id='chr',strand='+',start=1,end=903)
        with tempfile.TemporaryDirectory() as d:
            args=SimpleNamespace(flank_bp=100,genome='fake')
            with patch('orsearcher.engine.subprocess.run',return_value=SimpleNamespace(stdout='>region\n'+a['cds']+'\n')):
                values=initial_case(row,a,{},Path(d),args,{'chr':903},False)
            alt=next(v for v in values if v['source']=='discovery_exonerate')
            self.assertEqual(alt['alignment_tool'],'exonerate');self.assertEqual(alt['overlap'],1.)
    def test_unresolved_anchor_is_not_command_error(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);shard=root/'shard.fa';shard.write_text('>w\nATG\n');query=root/'queries.faa';query.write_text('>r\nMAAA\n')
            args=SimpleNamespace(translated_search='auto',jobs=1,threads=2,translated_evalue=.001,translated_timeout=10,
                translated_max_regions=5,workers=1,reference_metadata={},sample='test')
            shards=SimpleNamespace(shards=[shard],segments=[Segment('w','chr','',0,300,300,shard.name,False)],total_bases=300)
            def execute(cmd,output,**kwargs):Path(output).write_text('')
            class Evidence:
                def apply(self,models,label):pass
            with patch('orsearcher.translated.execute',execute),patch('orsearcher.translated.parse_hits',return_value=[h(1,90,1,30)]),patch('orsearcher.translated.reconstruct',return_value=([],dict(anchor='a',status='no_model',errors=[]))):
                self.assertEqual(discover(args,shards,{},[],root/'translated',Evidence(),query),[])
            import json
            summary=json.loads((root/'translated/summary.json').read_text())
            self.assertEqual(summary['command_failed_cases'],0);self.assertEqual(summary['unresolved_anchors'],1)
            self.assertIn('unresolved_discovery_anchors',summary['warnings'])
            self.assertNotIn('local_reconstruction_command_failures',summary['warnings'])
