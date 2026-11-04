import csv,shutil,subprocess,tempfile,unittest
from pathlib import Path
from orsearcher.fasta import FastaRecord,write_fasta


@unittest.skipUnless(shutil.which('Rscript'),'R/ape native tests require the server environment')
class NativeTreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        probe=subprocess.run(['Rscript','-e','quit(status=if(requireNamespace("ape",quietly=TRUE)) 0 else 1)'],
                             capture_output=True,text=True,timeout=30)
        if probe.returncode:
            raise unittest.SkipTest('Native tree tests require R package ape; run in the analysis environment')
    def check(self,gappy=False):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);tree='((o1:.1,o2:.1,o3:.1,cOR:.1):1,(g1:.1,g2:.1,g3:.1,cGR:.1):1);'
            (p/'tree.nwk').write_text(tree)
            values=[(k,'OR' if k.startswith('o') or k=='cOR' else 'GR','candidate' if k.startswith('c') else 'reference') for k in ['o1','o2','o3','g1','g2','g3','cOR','cGR']]
            (p/'ids.tsv').write_text('tree_id\toriginal_id\tfamily\tstate\n'+'\n'.join('\t'.join([k,k,f,s]) for k,f,s in values)+'\n')
            write_fasta(p/'aligned.faa',[FastaRecord(k,'',('-'*250+'A'*50) if gappy and s=='candidate' else 'M'+'A'*298+'C') for k,f,s in values])
            script=Path(__file__).resolve().parents[1]/'src/orsearcher/tree_qc.R'
            result=subprocess.run(['Rscript',str(script),str(p/'tree.nwk'),str(p/'ids.tsv'),str(p/'aligned.faa'),str(p)],capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stderr)
            with (p/'candidate_qc.tsv').open() as f:return {r['candidate_id']:r for r in csv.DictReader(f,delimiter='\t')}
    def test_good_reference_and_local_consensus_votes(self):
        rows=self.check();self.assertEqual(rows['cOR']['identity_tree_vote'],'OR');self.assertEqual(rows['cGR']['identity_tree_vote'],'GR')
    def test_gappy_alignment_does_not_vote_despite_tree_side(self):
        rows=self.check(True);self.assertEqual({r['identity_tree_vote'] for r in rows.values()},{'unresolved'})
