"""CLI-first, validated algorithm parameters. Defaults preserve rc3 decisions."""
import math

# name: (default, minimum, maximum, explanation); identity values are percent.
SPECS = {
 'adaptive_min_aa': (250,1,None,'Minimum adaptive search-seed protein length'),
 'adaptive_max_aa': (650,1,None,'Maximum adaptive search-seed protein length'),
 'adaptive_reference_coverage': (.7,0,1,'Minimum adaptive-seed reference coverage fraction'),
 'adaptive_min_identity': (20.,0,100,'Minimum adaptive-seed identity percent'),
 'first_seed_hits': (6,1,None,'Minimum Miniprot seed hits in the first discovery pass'),
 'second_seed_hits': (3,1,None,'Minimum Miniprot seed hits in the relaxed discovery pass'),
 'first_output_coverage': (.25,0,1,'First-pass Miniprot --outc coverage fraction'),
 'second_output_coverage': (.15,0,1,'Second-pass Miniprot --outc coverage fraction'),
 'first_alignment_ratio': (.3,0,1,'First-pass Miniprot -p secondary alignment ratio'),
 'second_alignment_ratio': (.2,0,1,'Second-pass Miniprot -p secondary alignment ratio'),
 'alignment_output_score': (.2,0,1,'Miniprot --outs relative output score'),
 'evidence_or_targets': (50,1,None,'Maximum retained external OR DIAMOND targets per protein'),
 'evidence_gr_targets': (5,1,None,'Maximum retained external GR DIAMOND targets per protein'),
 'evidence_evalue': (1e-3,0,10,'DIAMOND protein-evidence search E-value; must encompass identity thresholds'),
 'evidence_domain_overlap': (.25,0,1,'Minimum reference OR-domain overlap for non-curated homology evidence'),
 'candidate_min_aa': (80,1,None,'Minimum protein length admitted to structural reconstruction'),
 'candidate_profile_evalue': (.1,0,10,'Maximum OR-profile E-value for structural admission'),
 'candidate_bits': (35.,0,None,'Minimum homolog bitscore for structural admission'),
 'candidate_evalue': (1e-3,0,10,'Maximum homolog E-value for structural admission'),
 'candidate_query_coverage': (.35,0,1,'Minimum query coverage for structural admission'),
 'identity_evalue': (1e-5,0,10,'Maximum E-value for credible OR and GR homologs (same criteria)'),
 'identity_bits': (50.,0,None,'Minimum credible OR and GR homolog bitscore'),
 'identity_query_coverage': (.5,0,1,'Minimum credible OR and GR query coverage'),
 'identity_min_identity': (25.,0,100,'Minimum credible OR and GR amino-acid identity percent'),
 'distant_evalue': (1e-3,0,10,'Maximum E-value for distant OR and GR homologs'),
 'distant_bits': (35.,0,None,'Minimum distant OR and GR homolog bitscore'),
 'distant_query_coverage': (.35,0,1,'Minimum distant OR and GR query coverage'),
 'distant_min_identity': (15.,0,100,'Minimum distant OR and GR identity percent'),
 'profile_evalue': (1e-3,0,10,'Maximum profile domain E-value for OR/GR support'),
 'profile_coverage': (.25,0,1,'Minimum OR/GR profile coverage fraction'),
 'family_score_margin': (10.,0,None,'Required OR-minus-GR or GR-minus-OR bitscore margin'),
 'distant_score_margin': (15.,0,None,'Required OR-minus-GR margin for the distant-OR path'),
 'background_score_margin': (30.,0,None,'Required background-minus-best-family bitscore margin'),
 'complete_min_aa': (250,1,None,'Minimum protein length for complete coding models'),
 'complete_max_aa': (650,1,None,'Maximum protein length for complete coding models'),
 'partial_min_aa': (150,1,None,'Minimum accepted partial-model protein length'),
 'partial_max_aa': (650,1,None,'Maximum accepted partial-model protein length'),
 'model_hit_evalue': (1e-20,0,10,'Maximum homology E-value for final complete/partial model assessment'),
 'model_hit_bits': (100.,0,None,'Minimum homology bitscore for final model assessment'),
 'model_hit_identity': (30.,0,100,'Minimum homology identity percent for final model assessment'),
 'model_hit_coverage': (.7,0,1,'Minimum homology query coverage for final model assessment'),
 'partial_hit_identity': (35.,0,100,'Minimum homolog identity percent for accepted partial models'),
 'complete_hit_bits': (120.,0,None,'Minimum whole-homolog bitscore for completeness support'),
 'complete_query_coverage': (.8,0,1,'Minimum query coverage for whole-homolog completeness support'),
 'complete_reference_coverage': (.85,0,1,'Minimum subject coverage for whole-homolog completeness support'),
 'complete_terminal_slack_aa': (25.,0,None,'Minimum allowed terminal slack in whole-homolog alignment, aa'),
 'complete_terminal_slack_fraction': (.1,0,1,'Allowed terminal slack as fraction of reference length'),
 'repair_hit_bits': (80.,0,None,'Minimum external-reference bitscore to plan a polish repair'),
 'repair_hit_evalue': (1e-10,0,10,'Maximum external-reference E-value to plan a polish repair'),
 'repair_hit_identity': (30.,0,100,'Minimum external-reference identity percent to plan a polish repair'),
}
CHOICES = {
 'complete_evidence': ('qualified-template',('qualified-template','external-homology','coding-and-identity'),
  'Completeness evidence: qualified-template (rc3 default), external-homology (whole active external OR homolog), or coding-and-identity (intact coding plus OR identity; limited full-length confidence)'),
 'reconstruction_reference_policy': ('audited',('audited','active-external'),
  'Local reconstruction references: audited (default) or all active external OR proteins; self-generated search seeds never qualify'),
}
DEFAULTS={k:v[0] for k,v in SPECS.items()} | {k:v[0] for k,v in CHOICES.items()} | {'identity_use_tree':False,'adaptive_require_full_reference':True}

def reconstruction_metadata(args):
    data=getattr(args,'reference_metadata',{})
    if getattr(args,'reconstruction_reference_policy','audited')!='active-external':return data
    # This copy is for template SELECTION only. Never overwrite identity audit metadata.
    return {rid:dict(value,reconstruction_allowed=True) for rid,value in data.items() if not rid.startswith(('adaptive_','annotation_'))}

def add_arguments(parser):
 groups = {name:parser.add_argument_group(name) for name in (
     'Discovery and adaptive seed thresholds', 'Family identity and protein evidence thresholds',
     'Template selection and model assessment thresholds', 'Repair admission thresholds')}
 def group(name):
  if name.startswith(('first_','second_','alignment_','candidate_','adaptive_')):
   return groups['Discovery and adaptive seed thresholds']
  if name.startswith('repair_'):
   return groups['Repair admission thresholds']
  if name.startswith(('complete_','partial_','model_','reconstruction_')):
   return groups['Template selection and model assessment thresholds']
  return groups['Family identity and protein evidence thresholds']
 for name,(default,low,high,help_text) in SPECS.items():
  group(name).add_argument('--'+name.replace('_','-'),type=int if type(default) is int else float,default=default,
    help=help_text)
 for name,(default,choices,help_text) in CHOICES.items():
  group(name).add_argument('--'+name.replace('_','-'),choices=choices,default=default,help=help_text)
 import argparse
 group('adaptive_require_full_reference').add_argument('--adaptive-require-full-reference',action=argparse.BooleanOptionalAction,default=True,help='Require qualified full-reference support for internal adaptive search seeds')
 group('identity_use_tree').add_argument('--identity-use-tree',action=argparse.BooleanOptionalAction,default=False,
  help='Retain QC-approved diagnostic-tree votes during final identity assessment; default: disabled to preserve rc3 behavior. A tree alone does not establish OR identity.')

def validate(args):
 values={k:getattr(args,k,default) for k,default in DEFAULTS.items()}
 for k,(default,lo,hi,_) in SPECS.items():
  v=values[k]
  if not math.isfinite(v) or v<lo or (hi is not None and v>hi) or ('evalue' in k and v<=0):
   raise ValueError('Invalid --'+k.replace('_','-')+': '+str(v))
 for prefix in ('complete','partial','adaptive'):
  if values[prefix+'_min_aa']>values[prefix+'_max_aa']:raise ValueError(prefix+' minimum aa must not exceed maximum aa')
 if values['evidence_evalue']<max(values['identity_evalue'],values['distant_evalue'],values['model_hit_evalue'],values['candidate_evalue'],values['repair_hit_evalue']):
  raise ValueError('--evidence-evalue must encompass identity, distant, candidate, repair and final model E-values')
 args.thresholds=values
 return values

def policy(model=None,args=None):
 return DEFAULTS | (getattr(args,'thresholds',{}) if args is not None else (model or {}).get('_thresholds',{}))
