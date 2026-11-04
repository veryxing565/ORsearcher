import copy,csv,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from orsearcher.cli import parser
from orsearcher.tuning import validate,DEFAULTS
from orsearcher.workflow import dispatch,overlay,make_seeds
from orsearcher.polish import assessment
from orsearcher.identity import judge
from test_polish import model
from test_identity import protein,hit

class TuningTests(unittest.TestCase):
 def test_default_model_policy_is_unchanged(self):
  m=model();m['_thresholds']=dict(DEFAULTS);self.assertTrue(assessment(m,{'core'})['usable_complete'])
  m['homology_hits'][0]['ref']='other';self.assertFalse(assessment(m,{'core'})['usable_complete'])
 def test_external_homology_policy_still_requires_whole_alignment(self):
  m=model(reference='other');m['_thresholds']={'complete_evidence':'external-homology'}
  self.assertTrue(assessment(m,{'core'})['usable_complete'])
  m['homology_hits'][0]['se']=100;self.assertFalse(assessment(m,{'core'})['usable_complete'])
 def test_coding_policy_does_not_accept_damage_or_GR_conflict(self):
  m=model(reference='other');m['_thresholds']={'complete_evidence':'coding-and-identity'}
  self.assertTrue(assessment(m,{'core'})['usable_complete'])
  for field,value in [('has_stop',False),('canonical',False),('frameshifts',1),('fusion_pairs',['a'])]:
   damaged=copy.deepcopy(m);damaged[field]=value;self.assertFalse(assessment(damaged,{'core'})['usable_complete'])
  m['gr_homology_hits']=[dict(m['homology_hits'][0],bits=900,ref='gr')];m['gr_profile_domains']=[[{'domain_evalue':1e-30,'coverage':.9}]]
  self.assertFalse(assessment(m,{'core'})['usable_complete'])
 def test_identity_relaxation_is_symmetric_for_GR_competitors(self):
  m=protein();m['homology_hits']=[hit(bits=70,identity=20)];m['_thresholds']={'identity_min_identity':18}
  self.assertEqual(judge(m)['status'],'OR_supported')
  m['gr_homology_hits']=[hit('gr',bits=150,identity=20)];m['gr_profile_domains']=[[dict(domain_evalue=1e-10,coverage=.9)]]
  self.assertEqual(judge(m)['status'],'identity_unresolved')
 def test_invalid_numeric_parameters_rejected(self):
  for values in [['--identity-min-identity','101'],['--model-hit-evalue','nan'],['--complete-min-aa','700','--complete-max-aa','650'],['--evidence-evalue','1e-8']]:
   args=parser().parse_args(['run','--genome','g','--outdir','o']+values)
   with self.assertRaises(ValueError):validate(args)
 def test_help_exposes_registered_parameters(self):
  p=parser();help=p._subparsers._group_actions[0].choices['run'].format_help()
  for name in DEFAULTS:self.assertIn('--'+name.replace('_','-'),help)
 def test_overlay_replaces_nargs_and_boolean_flags(self):
  p=parser();tokens=overlay(['run','--genome','g','--outdir','o','--gr-hmm','a','b','--identity-use-tree'],['--gr-hmm','c','--no-identity-use-tree'],p)
  args=p.parse_args(tokens);self.assertEqual(args.gr_hmm,[Path('c')]);self.assertFalse(args.identity_use_tree)

class WorkflowTests(unittest.TestCase):
 def test_round_overrides_call_same_main_and_inherit(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp)/'run';calls=[]
   def normal_main(argv):
    args=parser().parse_args(argv);calls.append(args);args.outdir.mkdir(parents=True);return 0
   result=dispatch(['run','--genome','g','--outdir',str(root),'--rounds','3','--round-seeds','off','--identity-bits','40','--or-reference','first','--next-round','--or-reference','second','--identity-bits','35','--next-round','--candidate-min-aa','50'],parser(),normal_main)
   self.assertEqual(result,0);self.assertEqual([a.identity_bits for a in calls],[40,35,35]);self.assertEqual([a.or_reference for a in calls],[Path('first'),Path('second'),Path('second')]);self.assertEqual([a.candidate_min_aa for a in calls],[80,80,50]);self.assertTrue(all(a.rounds==1 for a in calls))
 def test_failed_round_stops_remaining_rounds(self):
  with tempfile.TemporaryDirectory() as tmp:
   calls=[]
   def fail(argv):calls.append(argv);return 1
   self.assertEqual(dispatch(['run','--genome','g','--outdir',str(Path(tmp)/'r'),'--rounds','3','--round-seeds','off'],parser(),fail),1)
   self.assertEqual(len(calls),1)
 def test_next_round_cannot_change_target_genome(self):
  with tempfile.TemporaryDirectory() as tmp:
   with self.assertRaises(ValueError):dispatch(['run','--genome','g','--outdir',str(Path(tmp)/'r'),'--rounds','2','--round-seeds','off','--next-round','--genome','other'],parser(),lambda argv:0)
 def test_auto_seeds_are_deduplicated_and_provenance_preserved(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=Path(tmp);results=root/'previous/samples/s/results';results.mkdir(parents=True)
   fields=['candidate_id','independent_locus','category','identity_status','completeness_blockers','has_start','has_stop','canonical_splicing','frameshifts','internal_stops','ambiguous_bases']
   with (results/'OR_catalogue.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
    for cid in ['a','b']:writer.writerow(dict(zip(fields,[cid,'True','complete','OR_supported','','True','True','True','0','0','0'])))
   (results/'OR_complete.proteins.faa').write_text('>a\nMAAAAA\n>b\nMAAAAA\n')
   out=root/'seeds.faa';self.assertEqual(make_seeds(root/'previous',out),1)
   provenance=json.loads(out.with_suffix('.provenance.json').read_text());self.assertEqual(len(provenance[0]['seed_loci']),2)

if __name__=='__main__':unittest.main()
