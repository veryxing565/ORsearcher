# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Bounded second-aligner checks for EVERY existing OR candidate.

Uses cached local genomic regions, never edits input models or fills sequence.
Does not use the old 'overlaps one other model => fusion' rule to select refs.
All alignments, including nearby ORs, remain available for the locus audit.
"""
import argparse
import concurrent.futures as futures
import json
import logging
import re
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from . import review_or_catalogue as review
from orsearcher.fasta import FastaRecord,read_fasta,write_fasta

LOG=logging.getLogger('whole_catalogue_crosscheck')


def miniprot_frames(path,offset,sequence):
    """Protein CIGAR units follow the miniprot 0.18 manual; validate span."""
    records=[]
    for line in path.read_text().splitlines():
        if not line.startswith('##PAF\t'):continue
        f=line.split('\t')[1:]
        tags={v.split(':',2)[0]:v.split(':',2)[2] for v in f[12:] if len(v.split(':',2))==3}
        cigar=tags.get('cg',''); strand=f[4]; pos=int(f[7]) if strand=='+' else int(f[8])
        original=pos; events=[]
        for size,op in re.findall(r'(\d+)([MIDFGNUV])',cigar):
            size=int(size)
            consumed=size*3 if op in 'MD' else 0 if op=='I' else size
            if op in 'FG':
                if strand=='+': coordinates=list(range(pos+offset+1,pos+offset+consumed+1))
                else: coordinates=list(range(pos+offset,pos+offset-consumed,-1))
                local=[p-offset-1 for p in coordinates]
                events.append(dict(coordinates=coordinates,operation=op,target_nt=consumed,
                    context_ambiguous=sum(b not in 'ACGT' for b in sequence[max(0,min(local)-100):max(local)+101])))
            pos+=consumed if strand=='+' else -consumed
        expected=int(f[8])-int(f[7])
        if abs(pos-original)!=expected:
            raise ValueError(f'Protein CIGAR target-span mismatch: {path}: {f[0]}')
        records.append(dict(reference=f[0],strand=strand,coverage=(int(f[3])-int(f[2]))/int(f[1]),
             start=int(f[7])+offset+1,end=int(f[8])+offset,events=events))
    return records


def vulgar_frames(path,offset,sequence):
    records=[]
    for line in path.read_text().splitlines():
        if not line.startswith('vulgar: '):continue
        f=line.split();strand=f[8];pos=int(f[6]);initial=pos;events=[]
        for i in range(10,len(f),3):
            op=f[i];q=int(f[i+1]);n=int(f[i+2])
            if op=='F':
                if strand=='+':coords=list(range(offset+pos+1,offset+pos+n+1))
                else:coords=list(range(offset+pos,offset+pos-n,-1))
                # A query-only frame gap has a boundary, not a fabricated base.
                boundary=offset+pos+(1 if strand=='+' else 0)
                loc=[p-offset-1 for p in coords] or [boundary-offset-1]
                events.append(dict(coordinates=coords,boundary=boundary,operation='F',query_aa=q,target_nt=n,
                   context_ambiguous=sum(b not in 'ACGT' for b in sequence[max(0,min(loc)-100):max(loc)+101])))
            pos+=n if strand=='+' else -n
        if pos!=int(f[7]):raise ValueError('Vulgar target-span mismatch')
        records.append(dict(reference=f[1],strand=strand,qstart=int(f[2])+1,qend=int(f[3]),
                      start=min(initial,pos)+offset+1,end=max(initial,pos)+offset,events=events))
    return records


def target_references(models,refs,previous,limit):
    chosen=next((m for m in models if m['id']==previous['selected_model_id']),models[0])
    ranking=sorted([m for m in models if m.get('ref') in refs and m.get('family')=='OR'],
          key=lambda m:(m['overlap']>=.5,m['ref_cov']>=.8,m['identity'],m['or_bits'],m['overlap']),reverse=True)
    names=[]
    # Keep the best whole-protein homolog in the structural test. Ranking only
    # the identity of a short local fit can silently omit the full-length hit.
    full=[r for r in chosen.get('full_ref_support',[]) if r in refs]
    if full and previous.get('model_category')!='pseudogene':names.append(full[0])
    # Highest-identity homolog tests disruption; a structurally intact alternate
    # needs its OWN homolog too, rather than only the same disruptive reference.
    if ranking and ranking[0]['ref'] not in names:names.append(ranking[0]['ref'])
    if full and full[0] not in names:names.append(full[0])
    alternatives=sorted([m for m in ranking if review.coding_intact(m) and m['canonical']],
                         key=lambda m:(m['has_start'] and m['has_stop'],m['identity'],m['ref_cov']),reverse=True)
    for m in alternatives:
        if m['ref'] not in names:names.append(m['ref']);break
    for m in ranking:
        if m['ref'] not in names:names.append(m['ref'])
    names.extend(r for r in refs if r not in names)
    distinct=[];seen=set()
    for name in names:
        sequence=refs[name].sequence
        if sequence in seen:continue
        distinct.append(name);seen.add(sequence)
        if len(distinct)==limit:break
    return chosen,distinct


def process(row,models,run,outdir,previous,args):
    cid=row['candidate_id'];case=outdir/'cases'/cid;case.mkdir(parents=True,exist_ok=True)
    cached=run/'cases'/cid
    seq=next(read_fasta(cached/'region.fna')).sequence
    meta=json.loads((cached/'region.json').read_text());offset=meta['offset']
    reference_path=cached/'references.faa'
    if not reference_path.exists() or not reference_path.stat().st_size:
        return [],dict(candidate_id=cid,status='no_reference',model_count=0)
    refs={r.id:r for r in read_fasta(reference_path)}
    overrides=getattr(args,'reference_overrides',None)
    if overrides is not None:refs=dict(overrides)
    limit=args.pseudogene_references if row['original_category']=='pseudogene' or previous['model_category']=='pseudogene' else args.references_per_candidate
    chosen,rids=target_references(models,refs,previous,limit)
    miniframes=miniprot_frames(cached/'miniprot.gff3',offset,seq)
    (case/'miniprot_frame_events.json').write_text(json.dumps(miniframes))
    original=next(m for m in models if m['source']=='original')
    # Include structurally conflicting alternative exons, not just the old model.
    bounds=[m for m in models if m['overlap']>=.5 and m.get('ref_cov',0)>=.7]
    bounds.append(chosen)
    left=max(offset+1,min(s for m in bounds for s,e in m['exons'])-5000)
    right=min(offset+len(seq),max(e for m in bounds for s,e in m['exons'])+5000)
    if right-left+1>args.max_window_bp:
        left=max(offset+1,min(s for s,e in chosen['exons'])-5000)
        right=min(offset+len(seq),max(e for s,e in chosen['exons'])+5000)
    if right-left+1>args.max_window_bp:return [],dict(candidate_id=cid,status='window_too_large')
    if not rids:return [],dict(candidate_id=cid,status='no_reference')
    targetseq=seq[left-offset-1:right-offset];new_offset=left-1
    write_fasta(case/'target.fna',[FastaRecord('region','',targetseq)])
    write_fasta(case/'queries.faa',[refs[r] for r in rids])
    (case/'region.json').write_text(json.dumps(dict(offset=new_offset,contig=row['sequence_id'],start=left,end=right)))
    cmd=[args.exonerate,'--model','protein2genome','--querytype','protein','--targettype','dna',
      '--query',str(case/'queries.faa'),'--target',str(case/'target.fna'),'--bestn','3','--score','100',
      '--proteinwordlen',str(getattr(args,'exonerate_word_length',3)),
      '--maxintron',str(args.max_intron),'--refine','region','--refineboundary','300','--dpmemory','128',
      '--showalignment','no','--showtargetgff','yes','--showvulgar','yes',
      '--ryo','CDS_BEGIN\n%tcs\nCDS_END\n']
    (case/'command.json').write_text(json.dumps(cmd));began=time.monotonic()
    fallback=False
    try:review.command(cmd,case/'exonerate.out',args.timeout)
    except subprocess.TimeoutExpired:
        result=dict(candidate_id=cid,status='timeout',seconds=round(time.monotonic()-began,2))
        (case/'status.json').write_text(json.dumps(result));return [],result
    except Exception as error:
        if 'exit=-11' not in str(error):
            result=dict(candidate_id=cid,status='error',error=str(error));(case/'status.json').write_text(json.dumps(result));return [],result
        retry=list(cmd);retry[retry.index('--refine')+1]='none'
        boundary=retry.index('--refineboundary');del retry[boundary:boundary+2]
        (case/'command.retry.json').write_text(json.dumps(retry))
        try:
            remaining=args.timeout-(time.monotonic()-began)
            if remaining<=0:raise subprocess.TimeoutExpired(retry,args.timeout)
            review.command(retry,case/'exonerate.out',remaining);fallback=True
        except Exception as retry_error:
            result=dict(candidate_id=cid,status='timeout' if isinstance(retry_error,subprocess.TimeoutExpired) else 'error',
                        error=str(retry_error),first_error=str(error),retry='without_refinement')
            (case/'status.json').write_text(json.dumps(result));return [],result
    frames=vulgar_frames(case/'exonerate.out',new_offset,targetseq)
    def parse_models(path,source):
        values=[]
        for a in review.parse_output(path):
            if a['strand']!=row['strand'] or not a['exons']:continue
            ref=refs[a['reference_id']]
            m=review.build_model(cid,source,targetseq,new_offset,a['exons'],row['strand'],
                    frameshifts=a['frameshifts'],ref=ref.id,qstart=a['qstart']+1,qend=a['qend'],
                    ref_length=len(ref.sequence),identity=a['identity'])
            if m['cds']!=a['reported_cds'].upper() and not m['frameshifts']:
                LOG.warning('Reject coordinate/CDS mismatch: %s %s',cid,ref.id);continue
            m.update(alignment_tool='exonerate',alignment=dict(qstart=a['qstart']+1,qend=a['qend'],qlen=len(ref.sequence)))
            m['overlap']=review.overlap_fraction(original['exons'],m['exons'])
            # Neighbours inform fusion checks, not automatic target replacement.
            values.append(m)
            for v in review.terminal_variants(m,targetseq,new_offset,refs):
                v['overlap']=review.overlap_fraction(original['exons'],v['exons']);values.append(v)
        return values
    extra=parse_models(case/'exonerate.out','exonerate')
    # A more sensitive seed may favour a different heuristic splice path.
    # Recover the complementary old path ONLY where external cross-tool
    # substantial structure does not agree. Both presets are still ONE tool.
    from .audit_or_loci import tool,external_alignment,structurally_agrees,endpoints_agree
    minis=[m for m in models if tool(m)=='miniprot' and external_alignment(m) and m['ref_cov']>=.8
           and m['canonical'] and review.coding_intact(m) and m['length']>=250 and m.get('family')=='OR']
    agreed=any(x['ref_cov']>=.8 and x['length']>=250 and x['canonical'] and review.coding_intact(x)
               and x['overlap']>=.5 and structurally_agrees(x,m) and endpoints_agree(x,m) for x in extra for m in minis)
    legacy_needed=getattr(args,'legacy_rescue','off')=='auto' and getattr(args,'exonerate_word_length',3)!=6 and not agreed
    legacy_status='not_needed';legacy_error='';legacy_fallback=False
    if legacy_needed and not getattr(args,'allow_legacy',True):legacy_status='deferred_budget'
    elif legacy_needed:
        legacy=list(cmd);legacy[legacy.index('--proteinwordlen')+1]='6'
        (case/'command.legacy.json').write_text(json.dumps(legacy))
        legacy_start=time.monotonic();legacy_limit=min(args.timeout,30)
        try:
            try:review.command(legacy,case/'exonerate_legacy.out',legacy_limit)
            except Exception as first_error:
                remaining=legacy_limit-(time.monotonic()-legacy_start)
                if 'exit=-11' not in str(first_error) or remaining<=0:raise
                retry=list(legacy);retry[retry.index('--refine')+1]='none'
                boundary=retry.index('--refineboundary');del retry[boundary:boundary+2]
                (case/'command.legacy.retry.json').write_text(json.dumps(retry))
                review.command(retry,case/'exonerate_legacy.out',remaining);legacy_fallback=True
            extra.extend(parse_models(case/'exonerate_legacy.out','exonerate_legacy'))
            frames.extend(vulgar_frames(case/'exonerate_legacy.out',new_offset,targetseq));legacy_status='completed'
        except Exception as error:
            legacy_status='timeout' if isinstance(error,subprocess.TimeoutExpired) else 'error';legacy_error=str(error)
    (case/'exonerate_frame_events.json').write_text(json.dumps(frames))
    extra=list({m['id']:m for m in extra}.values())
    (case/'models.json').write_text(json.dumps(extra))
    result=dict(candidate_id=cid,status='completed_with_warnings' if legacy_error else 'completed',seconds=round(time.monotonic()-began,2),
                references=len(rids),model_count=len(extra),region_bp=len(targetseq),fallback_without_refinement=fallback,
                external_structure_agreement=agreed,legacy_needed=legacy_needed,legacy_status=legacy_status,legacy_error=legacy_error)
    result['legacy_fallback_without_refinement']=legacy_fallback
    (case/'status.json').write_text(json.dumps(result));return extra,result


