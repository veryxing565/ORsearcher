import copy,json,tempfile,unittest,csv,signal,subprocess
from unittest.mock import patch,MagicMock
from pathlib import Path
from types import SimpleNamespace
from orsearcher.recall import meta_groups,parse_wise,robust_models,unique_anchors,export_candidates
from orsearcher.benchmark import evaluate
from orsearcher.cli import parser
from orsearcher.runtime import execute
from test_v21 import model

class RecallTests(unittest.TestCase):
    def meta(self,low=1,ref='R000001'):
        root=f'{ref}|seg|-|{low}'
        return '\n'.join([
            f'seg\tMetaEuk\tgene\t{low}\t900\t80\t-\t.\tTarget_ID={ref};TCS_ID={root}',
            f'seg\tMetaEuk\tmRNA\t{low}\t900\t80\t-\t.\tTCS_ID={root}_mRNA;Parent={root}',
            f'seg\tMetaEuk\texon\t601\t900\t40\t-\t.\tTCS_ID={root}_exon_0;Parent={root}_mRNA',
            f'seg\tMetaEuk\tCDS\t601\t900\t40\t-\t0\tTCS_ID={root}_CDS_0;Parent={root}_exon_0'])+'\n'
    def test_meta_hierarchy_not_exon_as_transcript(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.gff';p.write_text(self.meta());g=meta_groups(p)[0]
            self.assertEqual(g['exons'],[(601,900)]);self.assertEqual(g['strand'],'-');self.assertEqual(g['ref'],'R000001')
    def test_meta_suboptimal_sets_keep_low_coordinate_identifier(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.gff';p.write_text(self.meta(1)+self.meta(1000));self.assertEqual(len(meta_groups(p)),2)
    def test_meta_missing_parent_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.gff';p.write_text('seg\tMetaEuk\tCDS\t1\t30\t40\t+\t0\tTCS_ID=x;Parent=missing\n')
            with self.assertRaises(ValueError):meta_groups(p)
    def test_meta_nonzero_phase_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.gff';p.write_text(self.meta().replace('\t-\t0\t','\t-\t1\t'))
            with self.assertRaises(ValueError):meta_groups(p)
    def test_meta_native_unspecified_phase_requires_triplet_exons(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'m.gff';p.write_text(self.meta().replace('\t-\t0\t','\t-\t.\t'))
            self.assertEqual(meta_groups(p)[0]['exons'],[(601,900)])
            p.write_text(p.read_text().replace('\t601\t900\t','\t600\t900\t'))
            with self.assertRaises(ValueError):meta_groups(p)
    def test_wise_minus_coordinates_normalized(self):
        text='42.5 ref 4 100 region 900 601 1 0\nregion\tGeneWise\tmatch\t900\t601\t42\t-\t.\tg\nregion\tGeneWise\tcds\t900\t601\t0\t-\t0\tg\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'w';p.write_text(text);g=parse_wise(p)[0]
            self.assertEqual(g['exons'],[(601,900)]);self.assertEqual(g['indels'],1);self.assertTrue(g['query_span_qualified'])
    def test_wise_multigene_summary_is_not_each_gene_coverage(self):
        text='42 ref 1 100 region 1 900 0 1\nregion\tGeneWise\tmatch\t1\t300\t42\t+\t.\tg1\nregion\tGeneWise\tcds\t1\t300\t0\t+\t0\tg1\nregion\tGeneWise\tmatch\t601\t900\t42\t+\t.\tg2\nregion\tGeneWise\tcds\t601\t900\t0\t+\t0\tg2\n'
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'w';p.write_text(text);gs=parse_wise(p);self.assertEqual(len(gs),2);self.assertTrue(all(not g['query_span_qualified'] for g in gs))
    def test_wise_no_alignment_is_empty_not_pseudogene(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'w';p.write_text('Bits Query start end\n//\n');self.assertEqual(parse_wise(p),[])
    def test_fragment_is_not_a_robust_covered_mask(self):
        a=model();a['alignment_reference_qualified']=True
        b=copy.deepcopy(a);b['has_start']=False
        c=copy.deepcopy(a);c['alignment_reference_qualified']=False
        self.assertEqual(robust_models([a,b,c]),[a])
    def test_self_reference_cannot_mask_new_search_even_if_flagged(self):
        a=model();a.update(alignment_reference_qualified=True,ref='adaptive_0')
        b=copy.deepcopy(a);b['ref']='annotation_test'
        self.assertEqual(robust_models([a,b]),[])
    def test_anchor_reference_union_does_not_expand_or_bridge_spans(self):
        a=dict(contig='chr',strand='+',start=1,end=100,bits=100,eligible=True,references=['a'])
        b=dict(a,start=21,end=120,bits=90,references=['b'])
        c=dict(a,start=41,end=140,bits=80,references=['c'])
        result=unique_anchors([a,b,c])
        self.assertEqual(len(result),2);self.assertEqual(result[0]['end'],100)
        self.assertEqual(result[0]['references'],['a','b']);self.assertEqual(result[1]['references'],['c'])
        self.assertEqual(a['references'],['a'])
    def test_export_does_not_filter_GR_or_damage_or_mark_complete(self):
        a=model();a.update(contig='chr',family='GR_likely');b=copy.deepcopy(a);b.update(exons=[[1001,1903]],stops=[dict(position=20)])
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);s=export_candidates([a,b],out)
            with (out/'candidate_catalogue.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(s['candidate_model_records'],2);self.assertTrue(all(r['category']=='candidate' and r['confidence']=='unverified' for r in rows))
    def test_export_same_protein_different_loci_not_merged(self):
        a=model();a['contig']='chr';b=copy.deepcopy(a);b['exons']=[[1001,1903]]
        with tempfile.TemporaryDirectory() as d:self.assertEqual(export_candidates([a,b],Path(d))['candidate_model_records'],2)
    def test_exact_coordinate_duplicates_keep_tool_provenance(self):
        a=model();a['contig']='chr';b=copy.deepcopy(a);b['alignment_tool']='genewise'
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);s=export_candidates([a,b],out);self.assertEqual(s['candidate_model_records'],1)
            self.assertIn('genewise',(out/'candidate_catalogue.tsv').read_text())
    def test_candidate_benchmark_never_claims_complete(self):
        a=model();a['contig']='chr'
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);export_candidates([a],root/'pool')
            gff=root/'reference.gff';gff.write_text('chr\tref\tgene\t1\t903\t.\t+\t.\tID=g;Name=Or1a\nchr\tref\tmRNA\t1\t903\t.\t+\t.\tID=t;Parent=g\nchr\tref\tCDS\t1\t903\t.\t+\t0\tParent=t\n')
            s=evaluate(SimpleNamespace(annotation=gff,prediction=root/'pool',outdir=root/'benchmark'))
            self.assertEqual(s['OR_annotation_loci_recovered'],1);self.assertEqual(s['OR_loci_with_complete_splice_endpoint_match'],0)
            self.assertEqual(s['OR_loci_with_candidate_splice_endpoint_match'],1);self.assertTrue(s['records_are_not_gene_counts'])
    def test_cli_recall_budgets_available(self):
        a=parser().parse_args(['run','--genome','g','--or-reference','o','--gr-reference','r','--or-hmm','h','--gr-hmm','h2','--outdir','out','--recall-only'])
        self.assertTrue(a.recall_only);self.assertEqual(a.recall_search,'auto');self.assertEqual(a.recall_max_regions,500);self.assertEqual(a.recall_timeout,20)
    def test_meta_timeout_stops_own_process_group(self):
        with tempfile.TemporaryDirectory() as d:
            proc=MagicMock();proc.pid=12345;proc.wait.side_effect=[subprocess.TimeoutExpired('metaeuk',1),0,0]
            with patch('orsearcher.runtime.PROCESS_GROUP_SUPPORTED',True),patch('orsearcher.runtime.subprocess.Popen') as spawn,patch('orsearcher.runtime.os.killpg',create=True) as kill,patch('orsearcher.runtime.signal.SIGKILL',9,create=True):
                spawn.return_value.__enter__.return_value=proc
                output=Path(d)/'out'
                with self.assertRaises(RuntimeError):execute(['metaeuk','predictexons'],output,timeout=1,kill_process_group=True)
                self.assertTrue(spawn.call_args.kwargs['start_new_session'])
                self.assertEqual(kill.call_args_list[0].args,(12345,signal.SIGTERM))
                self.assertEqual(kill.call_args_list[1].args,(12345,signal.SIGKILL))
                status=json.loads(Path(str(output)+'.status.json').read_text())
                self.assertEqual(status['status'],'timeout');self.assertTrue(status['process_group_cleanup'])
