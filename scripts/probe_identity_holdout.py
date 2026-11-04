"""Exact-sequence-disjoint OR/GR/other protein assay; NOT gene annotation accuracy.

Curated labels and negative-background annotations are comparator labels, not
experimental truth. Held-out proteins cannot exactly match any reference DB.
Related homologs may remain; this is not a family/order holdout.
"""
import argparse,copy,importlib.util,json,random,time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from orsearcher.fasta import FastaRecord,read_fasta,write_fasta
from orsearcher.engine import prepare_references
from orsearcher.evidence import Evidence
from orsearcher.identity import judge,background_evidence
from orsearcher.runtime import save,file_digest
from orsearcher.support.review_or_catalogue import write_table


def select(records,count,rng):
    seen=set();items=[]
    for r in records:
        if r.sequence in seen or not 150<=len(r.sequence)<=1000:continue
        seen.add(r.sequence);items.append(r)
    items.sort(key=lambda r:r.id)
    if len(items)<count:raise ValueError('Insufficient distinct assay proteins')
    return rng.sample(items,count)


def panel(params,metadata,background,count,rng):
    return {family:select((r for r in read_fasta(Path(path)) if family!='OR' or metadata.get(r.id,{}).get('curated')),count,rng)
        for family,path in [('OR',params['or_reference']),('GR',params['gr_reference']),('other',background)]}


def run(args):
    begin=time.monotonic();args.outdir.mkdir(parents=True,exist_ok=False)
    params=json.loads(args.manifest.read_text())['parameters'];rng=random.Random(20261012)
    metadata=json.loads((args.manifest.parent/'shared/reference/reference_audit.json').read_text())
    groups=panel(params,metadata,args.background_reference,args.per_class,rng)
    queries=[FastaRecord(f'{family}_{i:03d}','original='+r.id,r.sequence)
        for family,items in groups.items() for i,r in enumerate(items)]
    query=args.outdir/'heldout.faa';write_fasta(query,queries);excluded={r.sequence for r in queries}
    background=args.outdir/'background_reference.faa'
    remaining=[r for r in read_fasta(args.background_reference) if r.sequence not in excluded]
    write_fasta(background,remaining)
    p=SimpleNamespace(**params)
    for field in ['or_reference','gr_reference','or_hmm']:setattr(p,field,Path(getattr(p,field)))
    p.gr_hmm=list(map(Path,p.gr_hmm));p.exclude_reference_pattern=[];p.exclude_reference_proteins=query
    p.threads=8;p.identity_policy='multichannel';p.background_reference=background
    refs=prepare_references(p,args.outdir);evidence=Evidence(p,args.outdir/'evidence',refs)
    # Confirm exact sequence removal from ALL three protein reference classes.
    for path in [Path(str(p.or_db)+'.faa'),Path(str(p.gr_db)+'.faa'),background]:
        assert not excluded & {r.sequence for r in read_fasta(path)},'Query/reference exact leakage'
    models=[dict(id=r.id,protein=r.sequence,length=len(r.sequence)) for r in queries]
    evidence.apply(models,'heldout_identity');background_evidence(models,p,args.outdir/'background')
    old=None
    if args.baseline_source:
        spec=importlib.util.spec_from_file_location('orsearcher.baseline_identity',args.baseline_source/'src/orsearcher/identity.py')
        old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    rows=[]
    for m,r in zip(models,queries):
        result=judge(m);family=r.id.split('_')[0]
        rows.append(dict(query_id=r.id,source_id=r.description,expected=family,status=result['status'],
            reason=result['reason'],baseline_status=old.judge(m)['status'] if old else '',**result['signals']))
    write_table(args.outdir/'identity_comparison.tsv',rows,list(rows[0]))
    result=dict(status='completed',seconds=round(time.monotonic()-begin,2),scope='protein_identity_sanity_not_genome_annotation_benchmark',
        exact_sequence_overlap_with_any_reference=0,per_class=args.per_class,
        counts={family:dict(Counter(r['status'] for r in rows if r['expected']==family)) for family in groups},
        baseline_counts={family:dict(Counter(r['baseline_status'] for r in rows if r['expected']==family)) for family in groups} if old else None,
        non_OR_assigned_OR=[r['query_id'] for r in rows if r['expected']!='OR' and r['status']=='OR_supported'],
        OR_assigned_non_OR=[r['query_id'] for r in rows if r['expected']=='OR' and r['status'].startswith('non_OR')],
        source_sha256={str(path):file_digest(path) for path in [args.manifest,args.background_reference,Path(params['or_reference']),Path(params['gr_reference'])]},
        limitations=['Small deterministic annotation-labelled protein panel, not an independent population test',
            'Related homologs remain; not family/order holdout','No genome structure or phylogenetic tree used','Unknown is not a correct negative family assignment'])
    save(args.outdir/'summary.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--background-reference',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--baseline-source',type=Path);p.add_argument('--per-class',type=int,default=32)
    run(p.parse_args())
