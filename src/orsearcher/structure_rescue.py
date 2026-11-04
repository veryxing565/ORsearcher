"""Bounded new-reference reconstruction, never splice or sequence synthesis."""
from .tuning import policy as tuning_policy, reconstruction_metadata
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from types import SimpleNamespace
from .identity import distant_hits,profile_signal,query_coverage
from .fasta import read_fasta
from .runtime import save,LOG
from .support import audit_or_loci as audit,review_or_catalogue as review,crosscheck_or_catalogue as cross


def reference_choices(models,refs,metadata,used_sequences,limit):
    scores={}
    for m in models:
        for h in distant_hits(m.get('homology_hits',[])):
            rid=h['ref'];info=metadata.get(rid,{})
            if (rid not in refs or not info.get('reconstruction_allowed')
                or refs[rid].sequence in used_sequences or not 100<=len(refs[rid].sequence)<=1500
                or h['bits']<35 or h['evalue']>1e-3 or query_coverage(h)<.35):continue
            cov=(h['se']-h['ss']+1)/max(1,h['slen'])
            score=(cov>=.6,cov>=.4,h.get('curated',False),h['bits'],cov,rid)
            scores[rid]=max(scores.get(rid,score),score)
    selected={};seen=set(used_sequences)
    for rid in sorted(scores,key=scores.get,reverse=True):
        seq=refs[rid].sequence
        if seq in seen:continue
        selected[rid]=refs[rid];seen.add(seq)
        if len(selected)>=limit:break
    return selected


def rescue_priority(models):
    """Only incomplete/ambiguous, credible OR-like coding paths are eligible."""
    original=next(m for m in models if m['source']=='original')
    local=[m for m in models if audit.near_target(m,original) and audit.structure_candidate(m)
           and not m.get('fusion_pairs') and not(m.get('repeated_or_domains') and m.get('repeated_reference_regions'))]
    if not local:return None
    decision=audit.decide_locus({},models,None)
    if decision['category']!='partial':return None
    hits=[h for m in local for h in distant_hits(m.get('homology_hits',[]))]
    if not hits:return None
    profile=any(profile_signal(m.get('or_profile_domains',[])) for m in local)
    high=any(h['bits']>=100 and h['evalue']<=1e-10 for h in hits)
    if not(profile or high):return None
    chosen=decision['model']
    cov=lambda h:(h['se']-h['ss']+1)/max(1,h['slen'])
    strongest=max((cov(h) for h in hits if query_coverage(h)>=.5),default=0)
    selected=max((cov(h) for h in distant_hits(chosen.get('homology_hits',[]))
                  if query_coverage(h)>=.5),default=0)
    span=max(e for s,e in original['exons'])-min(s for s,e in original['exons'])+1
    # Larger genomic windows with weak scattered fragments are late, rather
    # than using a raw short-query identity to win complementary compute.
    return (strongest-selected>=.1,strongest>=.5,profile,span<=100000,
            strongest,max(h['bits'] for h in hits),-span,original['cid'])


def legacy_priority(models):
    """Complementary seed budget goes to unresolved paths, not raw bit scores."""
    local=[m for m in models if audit.structure_candidate(m)]
    complete=any(review.complete_model(m) and audit.external_alignment(m) and m.get('ref_cov',0)>=.8 for m in local)
    disrupted=any(audit.external_alignment(m) and m.get('ref_cov',0)>=.9 and m.get('identity',0)>=.85
                  and (m.get('frameshifts') or m.get('stops')) for m in local)
    profile=any(profile_signal(m.get('or_profile_domains',[])) for m in local)
    return (disrupted,not complete,profile,
            max((m.get('ref_cov',0) for m in local if audit.external_alignment(m)),default=0),
            max((m.get('or_bits',0) for m in local),default=0))


def reconstruct(args,rows,byid,folder,evidence,refs):
    out=folder/'rescue';out.mkdir(exist_ok=True)
    budget=getattr(args,'structure_rescue_budget',0)
    if not budget:
        result=dict(status='disabled',selected=0,deferred=0,status_counts={})
        save(out/'summary.json',result);return result
    eligible=[];unavailable=[]
    for row in rows:
        cid=row['candidate_id'];priority=rescue_priority(byid[cid])
        if priority is None:continue
        cached=folder/'cases'/cid
        if not all((cached/p).exists() for p in ['region.fna','region.json','references.faa','miniprot.gff3']):
            unavailable.append(dict(candidate_id=cid,reason='missing_local_case'));continue
        old_query=folder/'crosscheck/cases'/cid/'queries.faa'
        used={r.sequence for r in read_fasta(old_query)} if old_query.exists() else set()
        chosen=reference_choices(byid[cid],refs,reconstruction_metadata(args),used,
                                 getattr(args,'structure_rescue_references',3))
        if not chosen:
            unavailable.append(dict(candidate_id=cid,reason='no_new_distinct_audited_reference'));continue
        eligible.append((priority,row,chosen))
    eligible.sort(key=lambda x:x[0],reverse=True)
    selected=eligible[:budget];deferred=eligible[budget:]
    save(out/'plan.json',dict(selected=[dict(candidate_id=r['candidate_id'],references=list(refs)) for _,r,refs in selected],
                             deferred=[r['candidate_id'] for _,r,_ in deferred],unavailable=unavailable,
                             target_annotation_used=False,not_a_new_discovery_search=True))
    statuses=[];extras=[]
    def job(row,chosen_refs):
        cid=row['candidate_id'];d=audit.decide_locus(row,byid[cid],None)
        params=SimpleNamespace(timeout=min(args.exonerate_timeout,getattr(args,'structure_rescue_timeout',20)),
            exonerate='exonerate',max_intron=args.max_intron,max_window_bp=args.max_window_bp,
            references_per_candidate=len(chosen_refs),pseudogene_references=len(chosen_refs),
            exonerate_word_length=6,legacy_rescue='off',allow_legacy=False,reference_overrides=chosen_refs)
        return cross.process(row,byid[cid],folder,out,
                             dict(selected_model_id=d['model']['id'],model_category=d['category']),params)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(job,row,chosen_refs):row for _,row,chosen_refs in selected}
        for future in as_completed(jobs):
            row=jobs[future]
            try:models,status=future.result()
            except Exception as error:models=[];status=dict(candidate_id=row['candidate_id'],status='error',error=str(error))
            for m in models:
                # Same true aligner, not an additional independent software vote.
                m.update(contig=row['sequence_id'],annotation_gene_ids=[],fusion_pairs=[],
                         source='exonerate_rescue'+('_terminal' if m['source'].endswith('_terminal') else ''),
                         id=row['candidate_id']+'.rescue_'+m['id'].rsplit('.',1)[-1])
            extras.extend(models);statuses.append(status)
    evidence.apply(extras,'bounded_alternative_reference_models')
    for m in extras:byid[m['cid']].append(m)
    result=dict(status='completed_with_warnings' if any(s['status']!='completed' for s in statuses) or deferred else 'completed',
                selected=len(selected),deferred=len(deferred),unavailable=len(unavailable),
                model_count=len(extras),status_counts=dict(Counter(s['status'] for s in statuses)),statuses=statuses,
                true_alignment_tool='exonerate',target_annotation_used=False)
    save(out/'summary.json',result);LOG.info('Bounded alternative-reference reconstruction: %s',
        {k:v for k,v in result.items() if k!='statuses'});return result
