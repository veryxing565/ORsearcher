import copy,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from orsearcher.identity import judge,structural_family
from orsearcher.refine import repair_rank,eligible,plan_groups
from orsearcher.engine import choose_refs
from orsearcher.cli import parser
from orsearcher.fasta import FastaRecord
from orsearcher.support.audit_or_loci import collect_neighbors
from test_v21 import model


def hit(ref='or1',bits=200,identity=50,coverage=1,curated=True):
    return dict(ref=ref,bits=bits,evalue=1e-30,identity=identity,qlen=300,slen=300,qs=1,qe=int(300*coverage),
        ss=1,se=int(300*coverage),curated=curated,reference_sequence_key=ref)


def protein():
    m=model();m['homology_hits']=[hit()];m['or_profile_domains']=[dict(domain_evalue=1e-20,coverage=.9)]
    m['gr_homology_hits']=[];m['gr_profile_domains']=[];return m


class IdentityTests(unittest.TestCase):
    def test_profile_alone_never_confirms_OR(self):
        m=protein();m['homology_hits']=[];self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_single_homolog_without_other_support_stays_unknown(self):
        m=protein();m['or_profile_domains']=[];self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_two_distinct_curated_homologs_can_support_profile_negative_OR(self):
        m=protein();m['or_profile_domains']=[];m['homology_hits']=[hit(),hit('or2')]
        self.assertEqual(judge(m)['status'],'OR_supported')
    def test_identical_reference_keys_not_two_curated_votes(self):
        m=protein();m['or_profile_domains']=[];h=hit('alias');h['reference_sequence_key']='or1';m['homology_hits'].append(h)
        self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_exact_qualified_reference_supports_identity_not_independence(self):
        m=protein();m['or_profile_domains']=[];m['qualified_reference_exact_matches']=['or1']
        self.assertEqual(judge(m)['status'],'OR_supported');self.assertNotIn('confidence',judge(m))
    def test_OR_profile_and_GR_competition_conflict_not_deleted(self):
        m=protein();m['gr_homology_hits']=[hit('gr1',250)];m['gr_profile_domains']=[[dict(domain_evalue=1e-20,coverage=.9)]]
        self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_supported_GR_excluded_from_OR_pool(self):
        m=protein();m['or_profile_domains']=[];m['gr_homology_hits']=[hit('gr1',250)]
        m['gr_profile_domains']=[[dict(domain_evalue=1e-20,coverage=.9)]]
        self.assertEqual(judge(m)['status'],'non_OR_GR_supported')
    def test_tree_side_alone_cannot_delete_OR(self):
        self.assertEqual(judge(protein(),{'tree_side':'GR_side','gap_fraction':.9})['status'],'OR_supported')
    def test_GR_tree_without_GR_competition_does_not_veto_OR(self):
        result=judge(protein(),{'identity_tree_vote':'GR'})
        self.assertEqual(result['status'],'OR_supported');self.assertTrue(result['signals']['tree_disagreement'])
    def test_reliable_OR_tree_can_support_weak_profile_negative_homology(self):
        m=protein();m['or_profile_domains']=[]
        self.assertEqual(judge(m,{'identity_tree_vote':'OR'})['status'],'OR_supported')
    def test_background_local_hit_is_not_a_deletion_rule(self):
        m=protein();m['or_profile_domains']=[];m['homology_hits']=[];m['background_hits']=[hit('neg',500,90,.3)]
        self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_near_full_reviewed_background_with_large_margin_can_exclude(self):
        m=protein();m['or_profile_domains']=[];m['homology_hits']=[];m['background_hits']=[hit('neg',500,90)]
        self.assertEqual(judge(m)['status'],'non_OR_other_supported')
    def test_strong_OR_profile_and_background_conflict_preserved(self):
        m=protein();m['background_hits']=[hit('neg',500,90)]
        self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_structural_family_preserves_old_family_audit(self):
        m=protein();m['family']='OR_GR_uncertain';structural_family(m)
        self.assertEqual(m['family'],'OR');self.assertEqual(m['legacy_family'],'OR_GR_uncertain')
    def test_weak_homolog_can_enter_repair_before_identity_confirmation(self):
        m=protein();m['or_evalue']=1;m['or_profile_domains']=[];m['homology_hits']=[hit(bits=45,coverage=.4)]
        self.assertTrue(eligible(m));self.assertEqual(judge(m)['status'],'identity_unresolved')
    def test_repair_rank_whole_external_support_not_longest_peptide(self):
        full=protein();long=copy.deepcopy(full);long['length']=900;long['full_ref_support']=[]
        self.assertGreater(repair_rank(full),repair_rank(long))
    def test_choose_refs_keeps_whole_homolog_and_distinct_sequences(self):
        m=protein();m['full_ref_support']=['whole'];m['homology_hits']=[hit('partial',300,.5)]
        m['search_references']=['alias','different'];refs={k:FastaRecord(k,'',seq) for k,seq in [('whole','M'*300),('partial','C'*300),('alias','M'*300),('different','D'*300)]}
        meta={k:dict(reconstruction_allowed=True) for k in refs}
        result=choose_refs(m,refs,3,meta);self.assertEqual(list(result),['whole','partial','different'])
    def test_missing_local_alignment_does_not_crash_neighbor_audit(self):
        with tempfile.TemporaryDirectory() as d:self.assertEqual(collect_neighbors([{'candidate_id':'c'}],Path(d),{}),[])
    def test_refine_CLI_and_budget_controls(self):
        args=parser().parse_args(['refine','--genome','g.fna','--or-reference','or.faa','--gr-reference','gr.faa','--or-hmm','o.hmm',
            '--gr-hmm','g.hmm','--outdir','new','--candidate-models','models.json','--background-reference','neg.faa'])
        self.assertEqual(args.command,'refine');self.assertEqual(args.structure_max_loci,500);self.assertEqual(args.annotation_mode,'unified')
    def test_structure_budget_only_counts_runnable_groups(self):
        args=SimpleNamespace(max_window_bp=500000,structure_references=3,structure_max_loci=1,reference_metadata={'r1':dict(reconstruction_allowed=True)})
        huge=protein();huge['exons']=[[1,300],[600000,600300]]
        absent=protein();absent.update(full_ref_support=[],homology_hits=[],ref='missing')
        good=protein();refs={'r1':FastaRecord('r1','','M'*300)}
        selected,deferred,unavailable=plan_groups([huge,absent,good],args,refs)
        self.assertEqual(selected,[good]);self.assertFalse(deferred)
        self.assertEqual([x['reason'] for x in unavailable],['window_too_large','no_qualified_reference'])
