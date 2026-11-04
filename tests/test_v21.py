import copy,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from orsearcher.references import supported_overlap
from orsearcher.evidence import route,Evidence
from orsearcher.engine import phylogeny
from orsearcher.fasta import FastaRecord,write_fasta
from orsearcher.support.audit_or_loci import decide_locus,tool
from orsearcher.support.crosscheck_or_catalogue import process
from orsearcher.support.finalize_catalogue_review import same_locus,export_unique
from orsearcher.benchmark import receptor_family,annotation_models


def model(source='original',ref='r1'):
    return dict(id='c.'+source+'.'+ref,cid='c',source=source,ref=ref,exons=[[1,903]],strand='+',
                cds='ATG'+'GCT'*299+'TAA',protein='M'+'A'*299,length=300,frameshifts=0,stops=[],
                ambiguous=0,remainder=0,has_start=True,has_stop=True,canonical=True,overlap=1.,bridge=False,
                family='OR',ref_cov=1.,identity=.9,or_ga=True,or_bits=200,or_evalue=1e-30,
                full_ref_support=['r1','r2'],reference_exact_matches=[],fusion_pairs=[],
                repeated_or_domains=False,repeated_reference_regions=False)


class V21Tests(unittest.TestCase):
    def test_tree_references_honor_prepared_exclusions(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            write_fasta(root/'or.faa',[FastaRecord('keptOR','','MAAA')])
            write_fasta(root/'gr.faa',[FastaRecord('keptGR','','MCCC')])
            write_fasta(root/'rawOR.faa',[FastaRecord('heldout_OR','','MAAA')])
            args=SimpleNamespace(or_db=root/'or',gr_db=root/'gr',or_reference=root/'rawOR.faa')
            self.assertEqual(phylogeny(args,[],{},{},root),{})
            text=(root/'phylogeny/id_map.tsv').read_text()
            self.assertIn('keptOR',text)
            self.assertNotIn('heldout_OR',text)

    def test_duplicate_full_references_are_one_sequence(self):
        e=Evidence.__new__(Evidence)
        protein='M'+'A'*299
        e.refs={r:FastaRecord(r,'',protein) for r in ['r1','r2']}
        e.sequence_references={protein:['r1','r2','unreported_exact']}
        e.reference_metadata={r:dict(curated=True,domains=[],reconstruction_allowed=True,sequence_key='same') for r in e.refs}
        e.hmm_ga=[(30,30),(30,30)]
        h=dict(identity=100,qlen=300,slen=300,qs=1,qe=300,ss=1,se=300,evalue=1e-50,bits=500)
        m=model()
        e.assign(m,dict(hmms=[[],[]],or_hits=[dict(h,ref=r) for r in e.refs],gr_hits=[]))
        self.assertEqual(len(m['full_ref_support']),1)
        self.assertEqual(m['curated_reference_count'],1)
        self.assertIn('unreported_exact',m['reference_exact_matches'])

    def test_benchmark_names_and_parent_mapping(self):
        self.assertEqual(receptor_family('Orco'),'OR')
        self.assertEqual(receptor_family('Or67d'),'OR')
        self.assertEqual(receptor_family('Gr21a'),'GR')
        self.assertEqual(receptor_family('odorant_binding_protein'),'')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'gold.gff3'
            p.write_text('chr\tRefSeq\tgene\t1\t903\t.\t+\t.\tID=g;gene=Or67d\n'
                         'chr\tRefSeq\tmRNA\t1\t903\t.\t+\t.\tID=t;Parent=g\n'
                         'chr\tRefSeq\tCDS\t1\t903\t.\t+\t0\tParent=t\n')
            genes,_,models,_=annotation_models(p)
            self.assertEqual(genes['g']['family'],'OR')
            self.assertEqual(models[0]['exons'],[(1,903)])

    def test_shared_splicing_not_proximity_merges_models(self):
        left=dict(exons=[[1,300],[401,700],[801,1100]])
        shared=dict(exons=[[201,300],[401,700],[801,1300]])
        unrelated=dict(exons=[[201,700],[801,1300]])
        self.assertTrue(same_locus(left,shared))
        self.assertFalse(same_locus(left,unrelated))

    def test_residual_overlap_downgrade_consistent_exports(self):
        with tempfile.TemporaryDirectory() as d:
            a=model();b=model('exonerate');b['exons']=[[302,1204]]
            rows=[dict(candidate_id=k,sequence_id='chr') for k in ['a','b']]
            decisions={k:dict(category='complete',confidence='reviewed_complete',model=m)
                       for k,m in [('a',a),('b',b)]}
            summary=export_unique(rows,decisions,Path(d))
            self.assertEqual(summary['residual_overlapping_pairs'],1)
            self.assertEqual(summary['confidence_counts'],{'provisional_complete':2})
            gff=(Path(d)/'nonredundant/OR_models.gff3').read_text()
            self.assertNotIn('confidence=reviewed_complete',gff)
            self.assertIn('locus_status=unresolved_overlap',gff)

    def test_non_OR_extension_is_not_domain_support(self):
        ds=[dict(start=227,end=428,domain_evalue=1e-30)]
        self.assertEqual(supported_overlap(dict(ss=1,se=214),ds),0)
        self.assertEqual(supported_overlap(dict(ss=227,se=428),ds),1)

    def test_homology_only_single_reference_not_reviewed(self):
        m=dict(or_bits=200,or_ga=False,gr_ga=False,or_evalue=1,family='OR_GR_uncertain',
               curated_reference_count=1,homology_hits=[])
        self.assertEqual(route(m),'low_priority')

    def test_original_is_miniprot_not_independent_tool(self):
        self.assertEqual(tool(model()),'miniprot')
        ms=[model(),model('miniprot','r1'),model('miniprot','r2')]
        self.assertEqual(decide_locus({},ms,None)['confidence'],'provisional_complete')
        ms.append(model('exonerate','r2'))
        self.assertEqual(decide_locus({},ms,None)['confidence'],'reviewed_complete')

    def test_exact_target_reference_not_independent(self):
        ms=[model(),model('miniprot','r1'),model('exonerate','r2')]
        for m in ms:m['reference_exact_matches']=['r1','r2']
        self.assertEqual(decide_locus({},ms,None)['confidence'],'provisional_complete')

    def test_original_and_miniprot_cannot_confirm_pseudogene(self):
        m=model();m['stops']=[dict(coordinates=[100,101,102])]
        damage=dict(kind='premature_stop',tools=['original','miniprot'],references=['r1','r2','r3'],
                    coordinates=[100,101,102],models=[m])
        self.assertEqual(decide_locus({},[m],damage)['category'],'partial')
        damage['tools']=['miniprot','exonerate']
        self.assertEqual(decide_locus({},[m],damage)['category'],'pseudogene')

    def test_longer_model_prevents_false_complete(self):
        ms=[model(),model('miniprot'),model('exonerate')]
        longer=model('miniprot','r3');longer.update(exons=[[1,1143]],length=380,ref_cov=.9,identity=.85)
        result=decide_locus({},ms+[longer],None)
        self.assertEqual(result['category'],'partial')
        self.assertIn('longer_same_locus',result['reason'])

    def fixture(self,root,empty=False):
        cached=root/'run/cases/c';cached.mkdir(parents=True)
        (cached/'region.fna').write_text('>region\n'+'ATG'+'GCT'*299+'TAA'+'\n')
        (cached/'region.json').write_text(json.dumps(dict(offset=0)))
        (cached/'references.faa').write_text('' if empty else '>r1\n'+'M'+'A'*299+'\n')
        return dict(candidate_id='c',sequence_id='chr',strand='+',original_category='complete')

    def test_empty_reference_returns_status(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=self.fixture(root,True)
            extra,status=process(row,[model()],root/'run',root/'out',{},SimpleNamespace())
            self.assertEqual(status['status'],'no_reference')

    def test_segfault_retries_without_refinement(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);row=self.fixture(root)
            args=SimpleNamespace(references_per_candidate=2,pseudogene_references=4,max_window_bp=500000,
                                 exonerate='exonerate',max_intron=200000,timeout=5)
            with patch('orsearcher.support.crosscheck_or_catalogue.miniprot_frames',return_value=[]), \
                 patch('orsearcher.support.crosscheck_or_catalogue.vulgar_frames',return_value=[]), \
                 patch('orsearcher.support.crosscheck_or_catalogue.review.parse_output',return_value=[]), \
                 patch('orsearcher.support.crosscheck_or_catalogue.review.command',side_effect=[RuntimeError('exit=-11'), '']) as cmd:
                _,status=process(row,[model()],root/'run',root/'out',dict(selected_model_id=model()['id'],model_category='complete'),args)
                self.assertTrue(status['fallback_without_refinement'])
                retry=cmd.call_args_list[1].args[0]
                self.assertEqual(retry[retry.index('--refine')+1],'none')


if __name__=='__main__':unittest.main()
