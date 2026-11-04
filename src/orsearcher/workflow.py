"""Repeat the normal CLI pipeline with inherited command-line overrides."""
import json,copy,hashlib
from pathlib import Path
from .fasta import read_fasta,write_fasta,FastaRecord
from .runtime import save

def option_groups(tokens,parser):
 actions=parser._subparsers._group_actions[0].choices['run']._option_string_actions
 result=[];i=0
 while i<len(tokens):
  token=tokens[i];name=token.split('=',1)[0]
  if name not in actions:
   if i==0 and token=='run':result.append(('command',[token]));i+=1;continue
   raise ValueError('Unknown round option: '+token)
  action=actions[name];group=[token];i+=1
  if '=' not in token and action.nargs!=0:
   if action.nargs in ('+','*'):
    while i<len(tokens) and tokens[i].split('=',1)[0] not in actions:group.append(tokens[i]);i+=1
   else:
    if i==len(tokens):raise ValueError('Missing value for '+name)
    group.append(tokens[i]);i+=1
  result.append((action.dest,group))
 return result

def overlay(base,overrides,parser):
 new=option_groups(overrides,parser);names={name for name,t in new}
 return [t for name,group in option_groups(base,parser) if name not in names for t in group]+[t for name,group in new for t in group]

def make_seeds(previous,out):
 import csv
 audit=[];records={}
 samples=sorted((previous/'samples').iterdir())
 for sample in samples:
  results=sample/'results'
  if not (results/'OR_catalogue.tsv').exists():raise ValueError('Previous round has no final catalogue: '+str(sample))
  with (results/'OR_catalogue.tsv').open() as handle:rows=list(csv.DictReader(handle,delimiter='\t'))
  proteins={}
  for kind in ['complete','partial']:
   f=results/f'OR_{kind}.proteins.faa'
   if f.exists() and f.stat().st_size:proteins.update({r.id:r.sequence for r in read_fasta(f)})
  selected=[]
  for r in rows:
   if r.get('independent_locus','True')!='True' or r['category'] not in ('complete','partial') or r.get('identity_status')!='OR_supported':continue
   if r.get('completeness_blockers','') not in ('','no_qualified_full_length_template_support'):continue
   if not all(r.get(k)=='True' for k in ['has_start','has_stop','canonical_splicing']):continue
   if not all(r.get(k)=='0' for k in ['frameshifts','internal_stops','ambiguous_bases']):continue
   sequence=proteins[r['candidate_id']];name='adaptive_round_'+hashlib.sha256(sequence.encode()).hexdigest()[:20]
   records[name]=FastaRecord(name,'self_generated search_only=true full_length_verified=false',sequence)
   selected.append(dict(seed_id=name,source_candidate=r['candidate_id'],category=r['category']))
  audit.append(dict(sample=sample.name,seed_loci=selected))
 # For a batch, do not let one sample's generated sequences seed another sample.
 if len(samples)!=1:raise ValueError('Automatic round seeds currently require one genome; use --round-seeds off for multi-sample manifests')
 write_fasta(out,records.values());save(out.with_suffix('.provenance.json'),audit)
 return len(records)

def dispatch(argv,parser,main):
 blocks=[[]]
 for token in argv:
  if token=='--next-round':blocks.append([])
  else:blocks[-1].append(token)
 first=parser.parse_args(blocks[0])
 if first.command!='run':raise ValueError('Multiple whole-pipeline rounds require the run command')
 if first.rounds<1:raise ValueError('--rounds must be positive')
 if len(blocks)>1 and len(blocks)!=first.rounds:raise ValueError('Provide exactly --rounds minus one --next-round blocks')
 if first.manifest and first.round_seeds!='off':raise ValueError('Multi-round manifests require --round-seeds off')
 for block in blocks[1:]:
  forbidden={'outdir','genome','manifest','sample','rounds','command','resume','help'}
  if any(name in forbidden for name,t in option_groups(block,parser)):raise ValueError('Next-round overrides cannot change target genome/sample/output')
 if first.resume:raise ValueError('Multi-round --resume is not supported yet; individual round directories support normal --resume')
 from .presets import set_run_defaults
 set_run_defaults(first)
 root=first.outdir.resolve()
 if root.exists() and any(root.iterdir()):raise FileExistsError('New multi-round output directory required: '+str(root))
 root.mkdir(parents=True,exist_ok=True)
 state=dict(status='running',rounds=first.rounds,completed=[],current=None,round_commands=[])
 save(root/'workflow_summary.json',state)
 inherited=overlay(blocks[0],['--rounds','1'],parser)
 try:
  for index in range(1,first.rounds+1):
   if index>1 and len(blocks)>1:
    forbidden={'outdir','genome','manifest','sample','rounds','command','resume','help'}
    if any(name in forbidden for name,t in option_groups(blocks[index-1],parser)):raise ValueError('Next-round overrides may change references, tuning and budgets, but not genome/sample/output identity')
    inherited=overlay(inherited,blocks[index-1],parser)
   folder=root/'rounds'/f'round{index:02d}'
   command=overlay(inherited,['--outdir',str(folder)],parser)
   parsed=parser.parse_args(command)
   if index>1 and parsed.round_seeds=='final-intact':
    seed=root/f'round{index:02d}_search_seeds.faa'
    count=make_seeds(root/'rounds'/f'round{index-1:02d}',seed)
    if count:command=overlay(command,['--search-seeds',str(seed)],parser)
    else:command=[t for name,group in option_groups(command,parser) if name!='search_seeds' for t in group]
   state['current']=index;state['round_commands'].append(command);save(root/'workflow_summary.json',state)
   code=main(command)
   if code:state.update(status='failed',error='Round '+str(index)+' failed');save(root/'workflow_summary.json',state);return code
   state['completed'].append(index);save(root/'workflow_summary.json',state)
  state.update(status='completed',current=None,final_output=str(folder));save(root/'workflow_summary.json',state);return 0
 except Exception as error:
  state.update(status='failed',error=str(error));save(root/'workflow_summary.json',state);raise
