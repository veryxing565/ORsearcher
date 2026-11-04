"""Annotation-blind, high-recall discovery; outputs are NOT confirmed OR genes.

MetaEuk raw exon sets are retained without its greedy overlap reduction.
GeneWise locally revisits weak/fragmentary anchors with audited external refs.
The target GFF is never opened here. Genomic bases are never corrected.
"""
from __future__ import annotations
from .tuning import policy as tuning_policy, reconstruction_metadata
import json,os,re,subprocess,time
from pathlib import Path
from collections import defaultdict,Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from types import SimpleNamespace
from .fasta import FastaRecord,read_fasta,write_fasta
from .models import parse_attributes
from .runtime import execute,save,digest,LOG
from .translated import parse_hits,anchors,chain_references,compact_models,reconstruct as complementary_reconstruct
from .support import review_or_catalogue as review
from .recall_schedule import exon_blocks,plan,region,reference_names,origin,rank


def meta_groups(path):
    """MetaEuk GFF has TCS_ID, and CDS -> exon -> mRNA -> gene parents."""
    features={};cds=[]
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t')
        if len(f)!=9:raise ValueError('Malformed MetaEuk GFF')
        a=parse_attributes(f[8]);key=a.get('ID',a.get('TCS_ID'))
        if not key:raise ValueError('MetaEuk feature lacks identifier')
        features[key]=(f,a)
        if f[2]=='CDS':cds.append(key)
    groups={}
    for key in cds:
        f,a=features[key];parent=key;seen=set()
        while features[parent][0][2]!='gene':
            if parent in seen:raise ValueError('Cyclic MetaEuk parents')
            seen.add(parent);parent=features[parent][1].get('Parent','')
            if parent not in features:raise ValueError('Missing MetaEuk parent')
        gf,ga=features[parent]
        if f[0]!=gf[0] or f[6]!=gf[6]:raise ValueError('Mixed MetaEuk contig/strand')
        s,e=int(f[3]),int(f[4])
        # Version 7 exports '.' for CDS phase. Its exon ORFs start/end at
        # codon boundaries; verify triplet length before treating '.' as zero.
        if s<1 or e<s or f[7] not in ['0','.'] or (e-s+1)%3:
            raise ValueError('Invalid MetaEuk CDS coordinate/phase')
        g=groups.setdefault(parent,dict(virtual_id=f[0],strand=f[6],ref=ga['Target_ID'],native_score=float(gf[5]),exons=[]))
        g['exons'].append((s,e))
    return list(groups.values())


def parse_wise(path):
    summaries=[];genes={}
    for line in Path(path).read_text().splitlines():
        f=line.split()
        if len(f)==9:
            try:
                summaries.append(dict(score=float(f[0]),ref=f[1],qs=int(f[2]),qe=int(f[3]),start=int(f[5]),end=int(f[6]),indels=int(f[7])))
                continue
            except ValueError:pass
        t=line.split('\t')
        if len(t)!=9 or t[1]!='GeneWise':continue
        if t[2]=='match':
            genes[t[8]]=dict(strand=t[6],exons=[],phase=0)
        elif t[2].lower() in ['cds','cds_exon']:
            if t[8] not in genes:raise ValueError('GeneWise CDS precedes match')
            g=genes[t[8]];s,e=sorted([int(t[3]),int(t[4])])
            if s<1:raise ValueError('Invalid GeneWise coordinate')
            if not g['exons']:g['phase']=int(t[7])
            g['exons'].append((s,e))
    if genes and len(summaries)!=1:raise ValueError('Ambiguous GeneWise query summary')
    out=[]
    for g in genes.values():
        if not g['exons']:continue
        out.append(dict(g,**summaries[0],query_span_qualified=len(genes)==1))
    return out


def robust_models(models):
    return [m for m in models if m.get('family')=='OR' and m.get('alignment_reference_qualified',False)
            and not m.get('ref','').startswith(('adaptive_','annotation_'))
            and m.get('full_ref_support') and m.get('ref_cov',0)>=.8 and m['length']>=250
            and m['has_start'] and m['has_stop'] and m['canonical'] and review.coding_intact(m)]


def unique_anchors(seeds):
    """Keep fixed coordinates, but retain refs from overlapping duplicate hints.

    Reference union must NOT expand the anchor span or transitively bridge loci.
    """
    ordered=sorted(seeds,key=lambda a:(0 if origin(a)=='translated_hsp' else 1,*rank(a)))
    unique=[];bins=defaultdict(list)
    for a in ordered:
        if not a['eligible']:continue
        ids={i for b in range(a['start']//100000,a['end']//100000+1) for i in bins[(a['contig'],a['strand'],b)]}
        match=next((i for i in sorted(ids) if max(0,min(a['end'],unique[i]['end'])-max(a['start'],unique[i]['start'])+1)
                    /max(a['end']-a['start']+1,unique[i]['end']-unique[i]['start']+1)>=.8),None)
        if match is not None:
            kept=unique[match]
            kept['references']=list(dict.fromkeys(kept['references']+a['references']))
            kept['origins']=sorted(set(kept.get('origins',[origin(kept)])+a.get('origins',[origin(a)])))
            continue
        i=len(unique);unique.append(dict(a,references=list(a['references']),origins=a.get('origins',[origin(a)])))
        for b in range(a['start']//100000,a['end']//100000+1):bins[(a['contig'],a['strand'],b)].append(i)
    return unique


def local_hints(args,raw,hsps,covered):
    """Build compact exon-block hints, keeping the original exon-set models."""
    seeds=anchors(hsps,robust_models(covered),args.reference_metadata)
    chain_references(seeds,hsps,args.max_intron,args.recall_window_bp-2*args.recall_flank_bp)
    for a in seeds:
        a.update(eligible=a['eligible'] or a['bits']>=30,origin='translated_hsp',
                 chain_coverage=max((r['query_coverage'] for r in a.get('reference_chains',[])),default=0))
    bins=defaultdict(list)
    for h in hsps:
        for b in range(h['start']//100000,h['end']//100000+1):bins[(h['contig'],h['strand'],b)].append(h)
    for m in raw:
        blocks=exon_blocks(m['exons'],args.recall_window_bp-2*args.recall_flank_bp)
        for i,block in enumerate(blocks):
            left=block[0][0];right=block[-1][1]
            nearby={}
            for b in range(left//100000,right//100000+1):
                for h in bins[(m['contig'],m['strand'],b)]:
                    if h['end']>=left and h['start']<=right:nearby[(h['ref'],h['start'],h['end'])]=h
            names=[h['ref'] for h in sorted(m.get('homology_hits',[]),key=lambda h:-h['bits'])]
            names.extend([m['ref']]);names.extend(h['ref'] for h in sorted(nearby.values(),key=lambda h:-h['bits']))
            seeds.append(dict(contig=m['contig'],strand=m['strand'],start=left,end=right,chain_start=left,chain_end=right,
                references=list(dict.fromkeys(names)),bits=m['native_score'],length=sum(e-s+1 for s,e in block)//3,
                eligible=True,origin='metaeuk',parent_model=m['id'],block_index=i,block_count=len(blocks),
                # A domain somewhere in a long parent is not evidence for
                # every extracted block. Parent scores are not broadcast.
                or_signal=len(blocks)==1 and (m.get('or_ga',False) or m.get('or_evalue',1)<=1e-3),
                or_bits=m.get('or_bits',0) if len(blocks)==1 else 0))
    return unique_anchors(seeds)


def oversized_hints(folder,mapping,refmap):
    """Retain local search opportunities from raw calls over the CDS limit.

    These light records are hints only, never exported as CDS/protein models.
    """
    values=[]
    for path in sorted(Path(folder).glob('*/skipped_models.json')):
        for item in json.loads(path.read_text()):
            g=item['model'];seg=mapping[g['virtual_id']];ref=refmap[g['ref']]
            exons=[[s+seg.original_start,e+seg.original_start] for s,e in g['exons']]
            values.append(dict(id='oversized_'+digest([seg.original_id,g['strand'],exons,ref])[:20],
                contig=seg.original_id,strand=g['strand'],ref=ref,exons=exons,native_score=g['native_score'],
                hint_only=True,raw_call_file=str(path),homology_hits=[]))
    return values


def rebuild(args,raw,hsps,covered,refs,lengths,folder,evidence):
    """Shared by fresh discovery and explicitly labelled cached-search tests."""
    unique=local_hints(args,raw,hsps,covered)
    selected,deferred,skipped,schedule=plan(args,unique,refs,lengths)
    save(folder/'anchors.json',unique);save(folder/'selected_anchors.json',selected)
    save(folder/'deferred_anchors.json',deferred);save(folder/'unschedulable_anchors.json',skipped);save(folder/'schedule_summary.json',schedule)
    LOG.info('Recall schedule hints=%d runnable=%d selected=%d deferred=%d skipped=%s',len(unique),schedule['runnable'],len(selected),len(deferred),schedule['skipped_counts'])
    values=[];statuses=[];by_index={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(wise_case,args,a,refs,lengths,folder/'genewise'):i for i,a in enumerate(selected)}
        for i,job in enumerate(as_completed(jobs),1):
            index=jobs[job]
            try:models,status=job.result()
            except Exception as error:models=[];status=dict(status='error',anchor=selected[index],error=str(error),errors=[str(error)])
            by_index[index]=models
            values.extend(models);statuses.append(status)
            if i%20==0 or i==len(selected):LOG.info('Recall local reconstruction %d/%d',i,len(selected))
    # Only genuinely unresolved cases consume complementary allowance. A full
    # looking alignment is a scheduling shortcut, not OR/fullness confirmation.
    needy=[i for i in range(len(selected)) if not any(substantial(m) for m in by_index[i])]
    budget=getattr(args,'recall_complementary_regions',100);chosen=needy[:budget]
    save(folder/'complementary_deferred.json',[selected[i] for i in needy[len(chosen):]])
    complementary=[]
    params=SimpleNamespace(**vars(args));params.translated_references=args.recall_references
    params.translated_flank_bp=args.recall_flank_bp;params.max_window_bp=args.recall_window_bp
    params.local_timeout=args.recall_timeout;params.exonerate_timeout=args.recall_timeout
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(complementary_reconstruct,params,selected[i],refs,folder/'complementary',lengths):i for i in chosen}
        for j,job in enumerate(as_completed(jobs),1):
            try:models,status=job.result()
            except Exception as error:models=[];status=dict(status='error',error=str(error),errors=[str(error)])
            for m in models:m['search_passes']=['recall_complementary']
            values.extend(models);complementary.append(dict(index=jobs[job],**status))
            if j%20==0 or j==len(chosen):LOG.info('Recall complementary reconstruction %d/%d',j,len(chosen))
    schedule['complementary']=dict(needed=len(needy),selected=len(chosen),deferred=len(needy)-len(chosen),statuses=complementary)
    save(folder/'schedule_summary.json',schedule)
    evidence.apply(values,'recall_rebuild_'+args.sample)
    save(folder/'new_local_models.json',values)
    return values,statuses,schedule


def substantial(m):
    return (m.get('ref_cov',0)>=.85 and m['length']>=250 and m['has_start'] and m['has_stop']
            and m['canonical'] and review.coding_intact(m))


def meta_shard(args,shard,mapping,query,refmap,refs,folder):
    case=folder/shard.stem;case.mkdir(parents=True,exist_ok=True)
    inputs=[shard,query];timeout=args.recall_search_timeout
    execute(['metaeuk','createdb',shard,case/'contigs','--dbtype','2'],case/'contigs.log',timeout=120,inputs=[shard])
    execute(['metaeuk','createdb',query,case/'references','--dbtype','1'],case/'references.log',timeout=120,inputs=[query])
    execute(['metaeuk','predictexons',case/'contigs',case/'references',case/'calls',case/'tmp',
             '--threads',max(1,args.threads//args.jobs),'-s',args.recall_sensitivity,'-e','100','--metaeuk-eval',str(args.recall_evalue),
             '--metaeuk-tcov','0.2','--min-length',args.recall_min_fragment_aa,'--min-exon-aa','11','--max-intron',args.max_intron,
             '--max-exon-sets',args.recall_exon_sets,'--max-seqs',max(1000,len(refmap)),
             '--split-memory-limit','16G','--mask','0','--comp-bias-corr','0'],case/'search.log',timeout=timeout,inputs=inputs,kill_process_group=True)
    execute(['metaeuk','unitesetstofasta',case/'contigs',case/'references',case/'calls',case/'predictions','--protein','1'],case/'export.log',timeout=300,inputs=inputs)
    if not Path(str(shard)+'.fai').exists():execute(['samtools','faidx',shard],case/'index.log',timeout=300,inputs=[shard])
    models=[];skipped=[]
    for g in meta_groups(case/'predictions.gff'):
        seg=mapping[g['virtual_id']];left=min(s for s,e in g['exons']);right=max(e for s,e in g['exons'])
        if right>seg.length:raise ValueError('MetaEuk CDS outside genome segment')
        if right-left+1>args.max_window_bp:
            skipped.append(dict(reason='window_too_large',model=g));continue
        result=subprocess.run(['samtools','faidx',str(shard),f"{g['virtual_id']}:{left}-{right}"],capture_output=True,text=True,check=True)
        seq=''.join(result.stdout.splitlines()[1:]).upper()
        if len(seq)!=right-left+1:raise ValueError('MetaEuk region length mismatch')
        rid=refmap[g['ref']];cid='meta_'+digest([seg.original_id,g['strand'],g['exons'],rid])[:16]
        m=review.build_model(cid,'metaeuk',seq,seg.original_start+left-1,[(s-left+1,e-left+1) for s,e in g['exons']],g['strand'],ref=rid,ref_length=len(refs[rid].sequence))
        # MetaEuk exon-set length is NOT measured reference coverage/identity.
        m.update(contig=seg.original_id,alignment_tool='metaeuk',ref_cov=0.,native_score=g['native_score'],protein_match=True,overlap=1.,
                 search_passes=['recall'],search_references=[rid],raw_alignment_count=1,annotation_gene_ids=[])
        models.append(m)
    save(case/'models.json',models);save(case/'skipped_models.json',skipped)
    return models,dict(shard=shard.name,status='completed',model_count=len(models),oversized_models=len(skipped))


def wise_case(args,a,refs,lengths,folder):
    key=digest([a['contig'],a['strand'],a['start'],a['end'],a['references']])[:16];case=folder/key;case.mkdir(parents=True,exist_ok=True)
    save(case/'anchor.json',a)
    def unresolved(reason):
        status=dict(anchor=key,status=reason,model_count=0,errors=[])
        save(case/'models.json',[]);save(case/'status.json',status);return [],status
    bounds=region(a,lengths[a['contig']],args.recall_flank_bp,args.recall_window_bp)
    if bounds is None:return unresolved('window_too_large')
    left,right=bounds
    p=subprocess.run(['samtools','faidx',str(args.genome),f"{a['contig']}:{left}-{right}"],capture_output=True,text=True,check=True)
    seq=''.join(p.stdout.splitlines()[1:]).upper()
    if len(seq)!=right-left+1:raise ValueError('GeneWise region length mismatch')
    write_fasta(case/'region.fna',[FastaRecord('region','',seq)])
    selected={rid:refs[rid] for rid in reference_names(a,refs,reconstruction_metadata(args),args.recall_references)}
    if not selected:return unresolved('no_reference')
    models=[];errors=[];begin=time.monotonic();env=os.environ.copy();env['WISECONFIGDIR']=str(args.wise_config)
    for i,(rid,ref) in enumerate(selected.items()):
        q=case/f'query_{i}.faa';output=case/f'wise_{i}.out';write_fasta(q,[ref])
        cmd=['genewise',q,case/'region.fna','-trev' if a['strand']=='-' else '-tfor','-alg','623','-init','local','-dymem','linear',
             '-sum','-gff','-quiet','-matrix','BLOSUM45.bla']
        save(case/f'command_{i}.json',list(map(str,cmd)))
        try:
            with output.open('w') as out,(case/f'wise_{i}.stderr.log').open('w') as err:
                p=subprocess.run(list(map(str,cmd)),stdout=out,stderr=err,env=env,timeout=args.recall_timeout)
            if p.returncode:raise RuntimeError(f'GeneWise exit {p.returncode}')
            for g in parse_wise(output):
                if g['strand']!=a['strand']:continue
                qs,qe=(g['qs'],g['qe']) if g['query_span_qualified'] else (0,0)
                m=review.build_model('wise_'+key,'genewise',seq,left-1,g['exons'],g['strand'],phase=g['phase'],frameshifts=g['indels'],
                                     ref=rid,qstart=qs,qend=qe,ref_length=len(ref.sequence),identity=0)
                if not g['query_span_qualified']:m['ref_cov']=0.
                if review.overlap_fraction([(a['start'],a['end'])],m['exons'])<.25:continue
                m.update(alignment_tool='genewise',alignment=dict(qstart=qs,qend=qe,qlen=len(ref.sequence)),native_score=g['score'])
                for v in [m,*review.terminal_variants(m,seq,left-1,selected)]:
                    v.update(contig=a['contig'],alignment_tool='genewise',protein_match=True,overlap=1.,search_passes=['recall'],
                             search_references=list(selected),raw_alignment_count=1,annotation_gene_ids=[])
                    models.append(v)
        except Exception as error:errors.append(dict(reference=rid,error=str(error)))
    save(case/'models.json',models)
    status=dict(anchor=key,status='completed_with_warnings' if errors else 'completed' if models else 'no_model',
                errors=errors,model_count=len(models),references=len(selected),seconds=round(time.monotonic()-begin,2))
    save(case/'status.json',status);return models,status


def discover(args,shards,refs,covered,folder,evidence,translated_folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    if args.recall_search=='off':save(folder/'summary.json',dict(status='disabled',warnings=[]));return []
    checkpoint=folder/'models.json'
    if checkpoint.exists() and (folder/'summary.json').exists() and json.loads((folder/'summary.json').read_text())['status']=='completed':return json.loads(checkpoint.read_text())
    start=time.monotonic();mapping={s.virtual_id:s for s in shards.segments};lengths={s.original_id:s.original_length for s in shards.segments}
    refmap={};records=[];seen=set()
    for rid,r in sorted(refs.items()):
        if r.sequence in seen:continue
        safe='R'+str(len(records)).zfill(6);refmap[safe]=rid;records.append(FastaRecord(safe,'',r.sequence));seen.add(r.sequence)
    query=folder/'references.faa';write_fasta(query,records);save(folder/'reference_map.json',refmap)
    raw=[];meta_status=[]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        jobs={pool.submit(meta_shard,args,s,mapping,query,refmap,refs,folder/'metaeuk'):s for s in shards.shards}
        for job in as_completed(jobs):
            try:values,status=job.result()
            except Exception as error:values=[];status=dict(shard=jobs[job].name,status='error',error=str(error))
            raw.extend(values);meta_status.append(status)
    evidence.apply(raw,'recall_meta_'+args.sample)
    save(folder/'metaeuk_models.json',raw)
    hsps=[]
    for path in Path(translated_folder).glob('*.hits.tsv'):hsps.extend(parse_hits(path,mapping))
    large=oversized_hints(folder/'metaeuk',mapping,refmap);save(folder/'oversized_local_hints.json',large)
    values,wise_status,schedule=rebuild(args,[*raw,*large],hsps,covered,refs,lengths,folder,evidence)
    raw.extend(values)
    evidence.apply(raw,'recall_all_'+args.sample)
    save(folder/'all_models.json',raw);models=compact_models(raw)
    warnings=[]
    if any(s['status']!='completed' for s in meta_status):warnings.append('metaeuk_failed_shards')
    if any(s.get('errors') for s in wise_status):warnings.append('genewise_failed_queries')
    if any(s.get('errors') for s in schedule['complementary']['statuses']):warnings.append('complementary_failed_queries')
    if schedule['complementary']['deferred']:warnings.append('complementary_budget_exceeded')
    if schedule['deferred']:warnings.append('genewise_budget_exceeded')
    if any(s['status']=='no_model' for s in wise_status):warnings.append('unresolved_genewise_anchors')
    if schedule['skipped_counts'].get('no_reference'):warnings.append('genewise_anchors_without_qualified_references')
    if schedule['skipped_counts'].get('window_too_large'):warnings.append('recall_unschedulable_windows')
    if any(s.get('oversized_models') for s in meta_status):warnings.append('metaeuk_oversized_models_retained_as_local_hints')
    save(folder/'summary.json',dict(status='completed_with_warnings' if warnings else 'completed',warnings=warnings,
         seconds=round(time.monotonic()-start,2),raw_models=len(raw),pooled_models=len(models),anchors=schedule['input_hints'],selected=schedule['selected'],deferred=schedule['deferred'],schedule=schedule,
         metaeuk_statuses=meta_status,genewise_statuses=wise_status,interpretation='Broad candidates, not confirmed OR identities or biological full-length genes'))
    save(checkpoint,models);return models


def export_candidates(models,folder):
    """No family/length/fullness filter, exact CDS dedup only; no gene counts."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);groups={};empty=0
    for m in models:
        if not m['exons'] or not m['protein']:empty+=1;continue
        key=(m['contig'],m['strand'],tuple(map(tuple,m['exons'])),m['cds'])
        g=groups.setdefault(key,dict(model=m,sources=set(),references=set()))
        g['sources'].add(m.get('alignment_tool',m['source']));g['references'].add(m.get('ref') or '')
    rows=[];proteins=[];cdss=[];gff=['##gff-version 3']
    for key,g in sorted(groups.items()):
        m=g['model'];cid='ORC_'+digest(key)[:20];exons=sorted(m['exons'],reverse=m['strand']=='-')
        row=dict(candidate_id=cid,sequence_id=m['contig'],strand=m['strand'],category='candidate',confidence='unverified',
                 protein_length=m['length'],sources=','.join(sorted(g['sources'])),references=','.join(sorted(g['references'])),
                 family_evidence=m.get('family','unknown'),or_bits=m.get('or_bits',0),gr_bits=m.get('gr_bits',0),
                 has_start=int(m['has_start']),has_stop=int(m['has_stop']),canonical=int(m['canonical']),internal_stops=len(m['stops']),
                 frameshift_hints=m['frameshifts'],interpretation='Unfiltered discovery model; not confirmed OR or confirmed pseudogene')
        rows.append(row);proteins.append(FastaRecord(cid,'unverified_discovery_candidate',m['protein']));cdss.append(FastaRecord(cid,'',m['cds']))
        left=min(s for s,e in exons);right=max(e for s,e in exons)
        gff.append('\t'.join(map(str,[m['contig'],'ORsearcher_recall','mRNA',left,right,'.',m['strand'],'.','ID='+cid+';category=candidate;confidence=unverified'])))
        n=0
        for i,(s,e) in enumerate(exons):
            gff.append('\t'.join(map(str,[m['contig'],'ORsearcher_recall','CDS',s,e,'.',m['strand'],(3-n%3)%3,f'ID={cid}.cds{i+1};Parent={cid}'])));n+=e-s+1
    review.write_table(folder/'candidate_catalogue.tsv',rows,list(rows[0]) if rows else ['candidate_id','sequence_id','strand','category','confidence'])
    write_fasta(folder/'candidate_proteins.faa',proteins);write_fasta(folder/'candidate_cds.fna',cdss);(folder/'candidate_models.gff3').write_text('\n'.join(gff)+'\n')
    summary=dict(candidate_model_records=len(rows),input_models=len(models),empty_models_not_exported=empty,exact_CDS_duplicates=len(models)-empty-len(rows),
                 family_evidence_counts=dict(Counter(r['family_evidence'] for r in rows)),not_gene_counts=True,not_confirmed_OR=True)
    save(folder/'candidate_summary.json',summary);return summary
