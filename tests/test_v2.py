import copy,json,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from orsearcher.discovery import union_models,adaptive_seeds
from orsearcher.evidence import route
from orsearcher.engine import choose_refs,annotation_references
from orsearcher.fasta import FastaRecord
from orsearcher.runtime import execute

class V2Tests(unittest.TestCase):
    def test_union_retains_first_only_and_separate_loci(self):
        m=dict(contig='c',strand='+',exons=[[1,9]],cds='ATGAAATAA',frameshifts=0,search_passes=['first'],search_references=['a'])
        second=dict(copy.deepcopy(m),search_passes=['second'],search_references=['b'])
        separate=dict(copy.deepcopy(second),contig='another')
        result=union_models([m],[second,separate])
        self.assertEqual(len(result),2);self.assertEqual(result[0]['search_passes'],['first','second'])
        first_only=union_models([copy.deepcopy(m)],[])
        self.assertEqual(len(first_only),1)
    def test_subga_not_discarded(self):
        m=dict(or_bits=80,homology_hits=[],or_evalue=1e-8,family='OR',gr_ga=False,or_ga=False)
        self.assertEqual(route(m),'OR_review')
        m.update(or_bits=0,or_evalue=1,family='GR_likely',gr_ga=True)
        self.assertEqual(route(m),'GR_likely')
        m.update(family='OR_GR_uncertain',gr_ga=False)
        self.assertEqual(route(m),'low_priority')
    def test_no_seed_no_annotation_is_valid(self):
        self.assertEqual(adaptive_seeds([]),[])
        self.assertEqual(annotation_references(SimpleNamespace(annotation=None),None,None),[])
    def test_reference_limit_and_lengths(self):
        # 2.2 budgets count distinct proteins, not aliases of one sequence;
        # long audited references may now enter repair, up to 1500 aa.
        refs={str(i):FastaRecord(str(i),'','M'+'ACDEFGHIKLMNPQR'[i]+'A'*298) for i in range(15)}
        refs['bad']=FastaRecord('bad','','M'*1600)
        result=choose_refs(dict(search_references=['bad',*refs]),refs,10)
        self.assertEqual(len(result),10);self.assertNotIn('bad',result)
    def test_failure_is_not_a_cache_hit(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'command.log';cmd=[sys.executable,'-c','raise SystemExit(3)']
            for _ in range(2):
                with self.assertRaises(RuntimeError):execute(cmd,out)
                self.assertEqual(json.loads(Path(str(out)+'.status.json').read_text())['status'],'error')
            self.assertFalse(out.exists())
    def test_success_cache_and_input_invalidation(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'command.log';inp=Path(d)/'input';inp.write_text('a')
            cmd=[sys.executable,'-c','print("ok")'];a=execute(cmd,out,inputs=[inp]);b=execute(cmd,out,inputs=[inp])
            self.assertEqual(a,b);inp.write_text('longer');c=execute(cmd,out,inputs=[inp])
            self.assertNotEqual(a['fingerprint'],c['fingerprint'])
