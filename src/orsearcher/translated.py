"""Independent genome-wide translated search followed by bounded reconstruction.

HSPs are discovery anchors, NOT CDS predictions or evidence of pseudogenization.
All coordinates are restored before pooling; neighboring/tandem anchors are
never joined by a transitive span merge. Target annotation is never read here.
"""
from __future__ import annotations
from .tuning import policy as tuning_policy, reconstruction_metadata
import concurrent.futures as futures
import json,subprocess
from collections import defaultdict,Counter
from pathlib import Path
from .fasta import FastaRecord,read_fasta,write_fasta
from .models import parse_miniprot_gff,attach_reference_lengths
from .sharding import Segment
from .runtime import execute,save,digest,LOG
from .support import review_or_catalogue as review

def parse_hits(path,mapping):
    hits=[]
    with Path(path).open() as f:
        for line in f:
            if not line.strip():continue
            q,s,identity,length,qs,qe,ss,se,evalue,bits,frame,qlen=line.rstrip().split('\t')
            seg=mapping[s];left,right=sorted([int(ss),int(se)])
            if left<1 or right>seg.length:raise ValueError('TBLASTN coordinate outside window')
            if int(length)<20 or float(bits)<30:continue
            hits.append(dict(contig=seg.original_id,strand='+' if int(frame)>0 else '-',start=left+seg.original_start,end=right+seg.original_start,
                ref=q,qstart=int(qs),qend=int(qe),qlen=int(qlen),identity=float(identity)/100,length=int(length),bits=float(bits),evalue=float(evalue)))
    # Exact boundary duplicates from overlapping genome windows add no evidence.
    return list({(h['contig'],h['strand'],h['start'],h['end'],h['ref'],h['qstart'],h['qend']):h for h in hits}.values())

def anchors(hits,covered,reference_metadata):
    index=defaultdict(list)
    for m in covered:
        for s,e in m['exons']:
            for b in range(s//100000,e//100000+1):index[(m['contig'],m['strand'],b)].append((s,e))
    groups=[];bins=defaultdict(list)
    for h in sorted(hits,key=lambda h:(h['bits'],h['length'],-h['evalue']),reverse=True):
        span=h['end']-h['start']+1
        if any(max(0,min(e,h['end'])-max(s,h['start'])+1)/span>=.8
               for b in range(h['start']//100000,h['end']//100000+1) for s,e in index[(h['contig'],h['strand'],b)]):continue
        candidates=set(i for b in range(h['start']//100000,h['end']//100000+1) for i in bins[(h['contig'],h['strand'],b)])
        match=None
        for i in sorted(candidates):
            a=groups[i][0]
            shared=max(0,min(a['end'],h['end'])-max(a['start'],h['start'])+1)
            if shared/min(a['end']-a['start']+1,span)>=.5:match=i;break
        if match is None:
            match=len(groups);groups.append([h])
            for b in range(h['start']//100000,h['end']//100000+1):bins[(h['contig'],h['strand'],b)].append(match)
        else:groups[match].append(h)
    out=[]
    for g in groups:
        h=dict(g[0]);names=list(dict.fromkeys(x['ref'] for x in g))
        keys={reference_metadata.get(n,{}).get('sequence_key',n) for n in names}
        # Weak single-reference matches are preserved but not expensively fitted.
        h.update(references=names,distinct_references=len(keys),eligible=h['bits']>=50 or h['length']>=25 and len(keys)>=2)
        out.append(h)
    return out

def chain_references(seeds,hits,max_intron,max_span):
    """Rank homologs by anchored, query-monotone HSP chains, not one exon.

    Chains are discovery hints only: genomic gaps may include unrelated DNA.
    A query-coordinate reset cannot join two tandem copies into a single chain.
    At most 160 strong HSPs/reference/window enter the quadratic local DP.
    """
    bins=defaultdict(list)
    for h in hits:
        for b in range(h['start']//100000,h['end']//100000+1):bins[(h['contig'],h['strand'],b)].append(h)
    for anchor in seeds:
        nearby={}
        radius=min(max_intron,max_span)
        for b in range(max(0,anchor['start']-radius)//100000,(anchor['end']+radius)//100000+1):
            for h in bins[(anchor['contig'],anchor['strand'],b)]:
                if h['end']<anchor['start']-radius or h['start']>anchor['end']+radius:continue
                nearby[(h['ref'],h['start'],h['end'],h['qstart'],h['qend'])]=h
        byref=defaultdict(list)
        for h in nearby.values():byref[h['ref']].append(h)
        ranked=[]
        for ref,hs in byref.items():
            hs=sorted(sorted(hs,key=lambda h:(-h['bits'],h['start'],h['qstart']))[:160],
                      key=lambda h:(h['start'] if h['strand']=='+' else -h['end'],h['qstart'],h['qend']))
            states=[]
            for i,h in enumerate(hs):
                g0,g1=(h['start'],h['end']) if h['strand']=='+' else (-h['end'],-h['start'])
                overlap=max(0,min(h['end'],anchor['end'])-max(h['start'],anchor['start'])+1)
                anchored=overlap/min(h['end']-h['start']+1,anchor['end']-anchor['start']+1)>=.5
                length=h['qend']-h['qstart']+1
                choices={anchored:dict(score=h['bits'],covered=length,path=[i],g0=g0,g1=g1)}
                for j,p in enumerate(hs[:i]):
                    pg1=p['end'] if p['strand']=='+' else -p['start']
                    gap=g0-pg1-1
                    if gap < -30 or gap>max_intron or h['qstart']<p['qend']-5 or h['qend']<=p['qend']:continue
                    novel=h['qend']-max(h['qstart']-1,p['qend'])
                    for touched,prev in states[j].items():
                        if g1-prev['g0']+1>max_span:continue
                        value=dict(score=prev['score']+h['bits']*novel/length-min(15,max(0,gap)/2000)-max(0,h['qstart']-p['qend']-1)*.08,
                                   covered=prev['covered']+novel,path=prev['path']+[i],g0=prev['g0'],g1=g1)
                        key=touched or anchored
                        if key not in choices or (value['score'],value['covered'])>(choices[key]['score'],choices[key]['covered']):choices[key]=value
                states.append(choices)
            anchored_states=[s[True] for s in states if True in s]
            if not anchored_states:continue
            best=max(anchored_states,key=lambda x:(x['score'],x['covered']))
            path=[hs[i] for i in best['path']]
            ranked.append(dict(ref=ref,score=round(best['score'],3),query_coverage=round(min(1,best['covered']/path[0]['qlen']),4),
                 hsp_count=len(path),start=min(h['start'] for h in path),end=max(h['end'] for h in path)))
        ranked.sort(key=lambda r:(-r['score'],-r['query_coverage'],r['ref']))
        anchor['reference_chains']=ranked
        # Hedge against a spurious multi-exon chain: two strong local references
        # still lead, then the coherent chains provide additional alternatives.
        anchor['references']=list(dict.fromkeys(anchor['references'][:2]+[r['ref'] for r in ranked]+anchor['references'][2:]))
        if any(r['hsp_count']>=2 and r['query_coverage']>=.25 and r['score']>=80 for r in ranked):anchor['eligible']=True
        best=next((r for r in ranked if r['hsp_count']>=2 and r['query_coverage']>=.25),None)
        anchor['chain_start']=best['start'] if best else anchor['start']
        anchor['chain_end']=best['end'] if best else anchor['end']
    return seeds

def reconstruct(args,anchor,refs,folder,lengths):
    key=digest([anchor['contig'],anchor['strand'],anchor['start'],anchor['end']])[:16]
    case=folder/'cases'/key;case.mkdir(parents=True,exist_ok=True)
    if (case/'models.json').exists() and (case/'status.json').exists() and json.loads((case/'status.json').read_text()).get('status')=='completed':
        return json.loads((case/'models.json').read_text()),json.loads((case/'status.json').read_text())
    # Recall scheduling may fall back from an overlong HSP chain to its
    # local anchor. Use the exact checked window for every reconstruction
    # tool, rather than re-expanding the rejected chain here.
    if 'planned_start' in anchor and 'planned_end' in anchor:
        left,right=anchor['planned_start'],anchor['planned_end']
    else:
        left=max(1,anchor.get('chain_start',anchor['start'])-args.translated_flank_bp);right=min(lengths[anchor['contig']],anchor.get('chain_end',anchor['end'])+args.translated_flank_bp)
    if not 1<=left<=right<=lengths[anchor['contig']]:raise ValueError('Invalid planned reconstruction window')
    save(case/'anchor.json',anchor)
    def unresolved(reason):
        status=dict(anchor=key,status=reason,model_count=0,errors=[])
        save(case/'models.json',[]);save(case/'status.json',status);return [],status
    if right-left+1>args.max_window_bp:return unresolved('window_too_large')
    p=subprocess.run(['samtools','faidx',str(args.genome),f"{anchor['contig']}:{left}-{right}"],capture_output=True,text=True,check=True)
    seq=''.join(p.stdout.splitlines()[1:]).upper();offset=left-1
    if len(seq)!=right-left+1:raise ValueError('Region length mismatch')
    query={n:refs[n] for n in anchor['references'] if n in refs and reconstruction_metadata(args).get(n,{}).get('reconstruction_allowed',False)}
    # Cap distinct sequences; chain-ranked alternatives compete with local hits.
    selected={};seen=set()
    for n,r in query.items():
        if r.sequence in seen:continue
        selected[n]=r;seen.add(r.sequence)
        if len(selected)==getattr(args,'translated_references',6):break
    if not selected:return unresolved('no_reconstruction_reference')
    write_fasta(case/'region.fna',[FastaRecord('region','',seq)]);write_fasta(case/'references.faa',selected.values())
    save(case/'region.json',dict(contig=anchor['contig'],offset=offset,end=right,strand=anchor['strand']))
    cid='translated_'+key;models=[];errors=[]
    cmd=['miniprot','-t','1','-k','5','-M','0','-L','15','-l','4','-n','2','--gff','-G',args.max_intron,'--outc','.15','--outs','.2','-N','12','--outn','6','-p','.2',case/'region.fna',case/'references.faa']
    try:
        execute(cmd,case/'miniprot.gff3',timeout=args.local_timeout,inputs=[case/'region.fna',case/'references.faa'])
        segment=Segment('region',anchor['contig'],'',offset,right,lengths[anchor['contig']],'region.fna',False)
        aligned=parse_miniprot_gff(case/'miniprot.gff3',{'region':segment});attach_reference_lengths(aligned,{n:len(r.sequence) for n,r in selected.items()})
        for c in aligned:
            if c.strand!=anchor['strand']:continue
            fs=sorted([f for f in c.features if f.feature_type=='CDS'],key=lambda f:f.start,reverse=c.strand=='-')
            if not fs:continue
            m=review.build_model(cid,'miniprot',seq,offset,[(f.start,f.end) for f in fs],c.strand,int(fs[0].phase) if fs[0].phase.isdigit() else 0,
                 c.frameshifts,c.reference_id,c.reference_start,c.reference_end,c.reference_length,c.identity)
            m.update(alignment_tool='miniprot',alignment=dict(qstart=c.reference_start,qend=c.reference_end,qlen=c.reference_length));models.append(m)
    except Exception as error:errors.append(dict(tool='miniprot',error=str(error)))
    near=lambda m:review.overlap_fraction([(anchor['start'],anchor['end'])],m['exons'])>=.5
    # Broad discovery is cheap; costly second reconstruction is only needed
    # when local Miniprot did not recover a plausible substantial coding model.
    substantial=any(near(m) and m['ref_cov']>=.8 and m['length']>=250 and review.coding_intact(m) and m['canonical']
                    and m['has_start'] and m['has_stop'] for m in models)
    if not substantial:
        cmd=['exonerate','--model','protein2genome','--querytype','protein','--targettype','dna','--query',case/'references.faa','--target',case/'region.fna',
             '--proteinwordlen',getattr(args,'exonerate_word_length',3),'--bestn','3','--score','60','--maxintron',args.max_intron,'--refine','none','--dpmemory','128','--showalignment','no','--showtargetgff','yes','--showvulgar','yes',
             '--ryo','CDS_BEGIN\n%tcs\nCDS_END\n']
        try:
            execute(cmd,case/'exonerate.out',timeout=args.exonerate_timeout,inputs=[case/'region.fna',case/'references.faa'])
            for a in review.parse_output(case/'exonerate.out'):
                if a['strand']!=anchor['strand'] or not a['exons']:continue
                r=selected[a['reference_id']]
                m=review.build_model(cid,'exonerate',seq,offset,a['exons'],a['strand'],frameshifts=a['frameshifts'],ref=r.id,qstart=a['qstart']+1,qend=a['qend'],ref_length=len(r.sequence),identity=a['identity'])
                if m['cds']!=a['reported_cds'].upper() and not m['frameshifts']:continue
                m.update(alignment_tool='exonerate',alignment=dict(qstart=a['qstart']+1,qend=a['qend'],qlen=len(r.sequence)));models.append(m)
        except Exception as error:errors.append(dict(tool='exonerate',error=str(error)))
    kept=[]
    for m in models:
        if not near(m):continue
        for v in [m,*review.terminal_variants(m,seq,offset,selected)]:
            v.update(contig=anchor['contig'],protein_match=True,overlap=1.,search_passes=['translated'],search_references=anchor['references'],
                     raw_alignment_count=1,annotation_gene_ids=[],alignment_tool=m['alignment_tool'])
            kept.append(v)
    status=dict(anchor=key,status='completed_with_warnings' if errors else 'completed' if kept else 'no_model',model_count=len(kept),miniprot_substantial=substantial,errors=errors,region_bp=len(seq),reference_count=len(selected))
    save(case/'models.json',kept);save(case/'status.json',status);return kept,status

def compact_models(models,ranking=None):
    """Quality-first nontransitive CDS anchors, retaining structural alternatives."""
    kept=[];bins=defaultdict(list)
    ordered=sorted(models,key=lambda m:tuple(map(str,(m['contig'],m['strand'],m['exons'],m.get('ref',''),m.get('alignment_tool',''),m['id']))))
    for m in sorted(ordered,key=ranking or review.model_rank,reverse=True):
        left=min(s for s,e in m['exons']);right=max(e for s,e in m['exons'])
        indices={i for b in range(left//100000,right//100000+1) for i in bins[(m['contig'],m['strand'],b)]}
        match=next((i for i in sorted(indices) if min(review.overlap_fraction(m['exons'],kept[i]['exons']),review.overlap_fraction(kept[i]['exons'],m['exons']))>=.7),None)
        if match is not None:
            parent=kept[match];alts=parent['discovery_alternatives']
            signature=lambda x:(tuple(map(tuple,x['exons'])),x.get('ref',''),x.get('alignment_tool',''),x['frameshifts'])
            if signature(m)!=signature(parent) and all(signature(m)!=signature(a) for a in alts):
                alts.append({k:v for k,v in m.items() if k!='discovery_alternatives'})
            continue
        i=len(kept);kept.append(dict(m,discovery_alternatives=[]))
        for b in range(left//100000,right//100000+1):bins[(m['contig'],m['strand'],b)].append(i)
    return kept

def discover(args,shards,refs,covered,folder,evidence,query):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    if args.translated_search=='off':
        save(folder/'summary.json',dict(status='disabled',warnings=[]));return []
    checkpoint=folder/'models.json'
    if checkpoint.exists() and (folder/'summary.json').exists():return json.loads(checkpoint.read_text())
    mapping={s.virtual_id:s for s in shards.segments}
    def search(shard):
        db=folder/(shard.stem+'.blast');marker=Path(str(db)+'.nhr')
        execute(['makeblastdb','-in',shard,'-dbtype','nucl','-parse_seqids','-out',db],folder/(shard.stem+'.db.log'),inputs=[shard],resume=marker.exists())
        output=folder/(shard.stem+'.hits.tsv')
        count=sum(s.shard==shard.name for s in shards.segments)
        execute(['tblastn','-query',query,'-db',db,'-num_threads',args.threads//args.jobs,'-evalue',args.translated_evalue,'-seg','yes','-soft_masking','true',
                 '-word_size','3','-max_target_seqs',max(1000,count),'-dbsize',shards.total_bases,'-max_intron_length','0',
                 '-outfmt','6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sframe qlen'],output,timeout=args.translated_timeout,inputs=[query,shard])
        return parse_hits(output,mapping)
    with futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:hsps=[h for block in pool.map(search,shards.shards) for h in block]
    hsps=list({(h['contig'],h['strand'],h['start'],h['end'],h['ref'],h['qstart'],h['qend']):h for h in hsps}.values())
    seed=anchors(hsps,covered,args.reference_metadata)
    chain_references(seed,hsps,getattr(args,'max_intron',200000),getattr(args,'max_window_bp',500000)-2*getattr(args,'translated_flank_bp',15000))
    eligible=[h for h in seed if h['eligible']]
    selected=eligible[:args.translated_max_regions];save(folder/'anchors.json',seed);save(folder/'deferred_anchors.json',eligible[len(selected):])
    LOG.info('Translated whole-genome search: HSPs=%d uncovered_anchors=%d eligible=%d selected=%d',len(hsps),len(seed),len(eligible),len(selected))
    lengths={s.original_id:s.original_length for s in shards.segments};raw=[];statuses=[]
    with futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(reconstruct,args,h,refs,folder,lengths):h for h in selected}
        for i,job in enumerate(futures.as_completed(jobs),1):
            try:models,status=job.result()
            except Exception as error:models=[];status=dict(anchor=jobs[job],status='error',error=str(error))
            raw.extend(models);statuses.append(status)
            if i%20==0 or i==len(selected):LOG.info('Translated local reconstruction %d/%d',i,len(selected))
    # Serial evidence cache access; worker threads never touch SQLite.
    evidence.apply(raw,'translated_'+args.sample)
    models=compact_models(raw);save(folder/'all_models.json',raw)
    warnings=[]
    if len(selected)<len(eligible):warnings.append('reconstruction_budget_exceeded')
    if any(s.get('errors') or s['status'] in ['error','completed_with_warnings'] for s in statuses):warnings.append('local_reconstruction_command_failures')
    unresolved=sum(s['status'] in ['no_model','no_reconstruction_reference','window_too_large'] for s in statuses)
    if unresolved:warnings.append('unresolved_discovery_anchors')
    save(folder/'summary.json',dict(status='completed_with_warnings' if warnings else 'completed',warnings=warnings,hsps=len(hsps),uncovered_anchors=len(seed),eligible=len(eligible),
         reconstructed=len(selected),deferred=len(eligible)-len(selected),raw_models=len(raw),pooled_models=len(models),
         unresolved_anchors=unresolved,command_failed_cases=sum(bool(s.get('errors')) or s['status']=='error' for s in statuses),
         preserved_alternatives=sum(len(m['discovery_alternatives']) for m in models),statuses=sorted(statuses,key=lambda s:str(s['anchor']))))
    save(checkpoint,models)
    return models
