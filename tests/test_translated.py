import tempfile,unittest,copy
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace
from orsearcher.translated import parse_hits,anchors,compact_models,discover
from orsearcher.sharding import Segment
from orsearcher.support.audit_or_loci import tool
from orsearcher.support.review_or_catalogue import terminal_variants
from test_v21 import model

def hit(start,end,ref='r1',bits=80):
    return dict(contig='chr',strand='+',start=start,end=end,ref=ref,bits=bits,length=(end-start+1)//3,evalue=1e-8)

class TranslatedTests(unittest.TestCase):
    def test_remap_strand_and_boundary_duplicates(self):
        mapping={'w1':Segment('w1','chr','',100,1000,2000,'s',True),'w2':Segment('w2','chr','',200,1100,2000,'s',True)}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits';p.write_text('r1\tw1\t30\t30\t1\t30\t201\t290\t1e-8\t60\t1\t400\n'
                'r1\tw2\t30\t30\t1\t30\t101\t190\t1e-8\t60\t1\t400\n'
                'r2\tw1\t30\t30\t1\t30\t490\t401\t1e-8\t60\t-2\t400\n')
            hs=sorted(parse_hits(p,mapping),key=lambda h:h['start'])
            self.assertEqual(len(hs),2);self.assertEqual((hs[0]['start'],hs[0]['end'],hs[0]['strand']),(301,390,'+'))
            self.assertEqual((hs[1]['start'],hs[1]['end'],hs[1]['strand']),(501,590,'-'))
    def test_invalid_coordinate_fails(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'hits';p.write_text('r1\tw\t30\t30\t1\t30\t1\t290\t1e-8\t60\t1\t400\n')
            with self.assertRaises(ValueError):parse_hits(p,{'w':Segment('w','chr','',0,100,100,'s',False)})
    def test_tandem_anchors_not_joined_transitively(self):
        hs=[hit(1,90),hit(151,240),hit(80,160,bits=40)]
        self.assertEqual(len(anchors(hs,[],{})),3)
    def test_covered_coding_not_whole_span_is_skipped(self):
        covered=[dict(contig='chr',strand='+',exons=[[1,90],[501,590]])]
        hs=anchors([hit(1,90),hit(201,290)],covered,{})
        self.assertEqual(len(hs),1);self.assertEqual(hs[0]['start'],201)
    def test_duplicate_reference_ids_not_independent(self):
        hs=anchors([hit(1,90,'r1',40),hit(1,90,'r2',40)],[],{'r1':{'sequence_key':'same'},'r2':{'sequence_key':'same'}})
        self.assertEqual(hs[0]['distinct_references'],1);self.assertFalse(hs[0]['eligible'])
    def test_true_aligner_survives_original_slot(self):
        self.assertEqual(tool(dict(source='original',alignment_tool='exonerate')),'exonerate')
        m=model();m['alignment_tool']='exonerate';m['cds']='ATG'+'GCT'*299;m['exons']=[[1,900]];m['has_stop']=False
        variants=terminal_variants(m,m['cds']+'TAA',0,{})
        self.assertTrue(variants);self.assertTrue(all(tool(v)=='exonerate' for v in variants))
    def test_compact_keeps_distinct_tandem_loci(self):
        a=model();a['contig']='chr';b=copy.deepcopy(a);b['exons']=[[1000,1902]]
        self.assertEqual(len(compact_models([a,copy.deepcopy(a),b])),2)
    def test_disabled_search_explicit_status(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(discover(SimpleNamespace(translated_search='off'),None,None,None,Path(d),None,None),[])
            self.assertIn('disabled',(Path(d)/'summary.json').read_text())
    def test_empty_live_search_contract_and_supported_blast_option(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);shard=root/'shard.fa';shard.write_text('>w\nATG\n');query=root/'queries.faa';query.write_text('>r\nMAAA\n')
            args=SimpleNamespace(translated_search='auto',jobs=1,threads=2,translated_evalue=.001,translated_timeout=10,
                translated_max_regions=5,workers=1,reference_metadata={},sample='test')
            shards=SimpleNamespace(shards=[shard],segments=[Segment('w','chr','',0,3,3,shard.name,False)],total_bases=3)
            calls=[]
            def execute(cmd,output,**kwargs):
                calls.append(cmd);Path(output).write_text('');return {'status':'completed'}
            class Evidence:
                def apply(self,models,label):self.models=models
            evidence=Evidence()
            with patch('orsearcher.translated.execute',execute):
                self.assertEqual(discover(args,shards,{},[],root/'translated',evidence,query),[])
            blast=next(c for c in calls if c[0]=='tblastn')
            self.assertEqual(blast[blast.index('-max_intron_length')+1],'0')
            self.assertEqual(blast[blast.index('-dbsize')+1],3)
            self.assertEqual(blast[blast.index('-max_target_seqs')+1],1000)
            self.assertEqual(evidence.models,[])
