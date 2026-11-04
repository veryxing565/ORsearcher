import copy,unittest
from test_identity import hit
from test_v21 import model
from orsearcher.identity import judge
from orsearcher.support.audit_or_loci import partial_selection_evidence,partial_selection_rank,decide_locus,structural_homology_hits


def local_hsp(query_length,start=1,bits=62.8):
    h=hit('shared',bits=bits,identity=24.6,coverage=1)
    h.update(qlen=query_length,slen=378,qs=start,qe=start+129,ss=227,se=356,evalue=1e-12)
    return h


def partial(source='original',length=345,start=True,stop=True):
    m=model(source);m.update(length=length,has_start=start,has_stop=stop,full_ref_support=[],or_bits=62.8,
        id=source,homology_hits=[local_hsp(length,209 if length==345 else 41)],or_profile_domains=[],
        gr_homology_hits=[],gr_profile_domains=[])
    return m


class Candidate10Tests(unittest.TestCase):
    def test_same_hsp_not_erased_on_longer_query(self):
        long=partial();short=partial('exonerate',176,False,False)
        a=partial_selection_evidence(long,[long,short]);b=partial_selection_evidence(short,[long,short])
        self.assertEqual(a['reference_coverage'],b['reference_coverage'])
        self.assertEqual(a['shared_bits'],b['shared_bits'])
        self.assertLess(a['query_fraction_diagnostic'],.5);self.assertGreater(b['query_fraction_diagnostic'],.5)
        self.assertFalse(a['query_fraction_used_in_rank'])
        self.assertGreater(partial_selection_rank(long,[long,short]),partial_selection_rank(short,[long,short]))
        self.assertIs(decide_locus({},[long,short],None)['model'],long)

    def test_below_old_distant_query_fraction_still_comparable(self):
        long=partial(length=600);long['homology_hits']=[local_hsp(600,350)]
        short=partial('exonerate',176,False,False)
        self.assertEqual(partial_selection_evidence(long,[long,short])['shared_bits'],62.8)

    def test_query_length_alone_does_not_change_rank(self):
        a=partial();b=copy.deepcopy(a);b['length']=600;b['homology_hits'][0]['qlen']=600
        self.assertEqual(partial_selection_rank(a,[a,b]),partial_selection_rank(b,[a,b]))

    def test_reference_aliases_compare_by_sequence_key(self):
        a=partial();b=copy.deepcopy(a);b['homology_hits'][0]['ref']='alias'
        self.assertEqual(partial_selection_rank(a,[a,b]),partial_selection_rank(b,[a,b]))

    def test_query_position_shift_with_same_reference_span_no_score_effect(self):
        a=partial();b=copy.deepcopy(a);b['homology_hits'][0].update(qs=20,qe=149)
        self.assertEqual(partial_selection_rank(a,[a,b]),partial_selection_rank(b,[a,b]))

    def test_native_identity_from_another_reference_not_tie_breaker(self):
        a=partial();b=copy.deepcopy(a);b.update(identity=.99,ref='different',ref_cov=.99)
        self.assertEqual(partial_selection_rank(a,[a,b]),partial_selection_rank(b,[a,b]))

    def test_better_shared_reference_support_can_select_shorter(self):
        a=partial();b=partial('exonerate',176,False,False)
        b['homology_hits'][0].update(bits=110,ss=170,se=356)
        self.assertGreater(partial_selection_rank(b,[a,b]),partial_selection_rank(a,[a,b]))

    def test_original_without_shared_support_not_protected(self):
        a=partial();a['homology_hits']=[]
        b=partial('exonerate',176,False,False)
        self.assertGreater(partial_selection_rank(b,[a,b]),partial_selection_rank(a,[a,b]))

    def test_strong_single_local_hit_not_new_identity_confirmation(self):
        m=partial(length=600);m['homology_hits']=[local_hsp(600,300,bits=300)]
        self.assertTrue(structural_homology_hits(m));self.assertEqual(judge(m)['status'],'identity_unresolved')

    def test_insignificant_hit_does_not_become_structural_support(self):
        m=partial();m['homology_hits'][0]['evalue']=.1
        self.assertFalse(structural_homology_hits(m))

    def test_tiny_reference_match_not_structural_support(self):
        m=partial();m['homology_hits'][0].update(se=250)
        self.assertFalse(structural_homology_hits(m))

    def test_damaged_long_anchor_not_forced_over_intact_short_path(self):
        a=partial();a['frameshifts']=1
        b=partial('exonerate',176,False,False)
        self.assertIs(decide_locus({},[a,b],None)['model'],b)

    def test_noncanonical_long_anchor_not_forced_over_canonical_path(self):
        a=partial();a['canonical']=False
        b=partial('exonerate',176,False,False)
        self.assertIs(decide_locus({},[a,b],None)['model'],b)

    def test_fused_long_anchor_not_restored_by_preservation(self):
        a=partial();a['fusion_pairs']=[['left','right']]
        b=partial('exonerate',176,False,False)
        self.assertIs(decide_locus({},[a,b],None)['model'],b)

    def test_endpoint_tie_break_does_not_upgrade_partial_to_complete(self):
        a=partial();b=partial('exonerate',176,False,False)
        result=decide_locus({},[a,b],None)
        self.assertEqual(result['category'],'partial');self.assertEqual(result['confidence'],'unresolved')
