"""Structure first, family identity second; preserve the broad input pool.

Uses no target annotation to decide inclusion or select a model. Repair budget
is explicit. Unknown identity is NOT converted into evidence of a non-OR.
"""
from pathlib import Path
from types import SimpleNamespace
from collections import Counter
import copy,json,time
from .runtime import save,digest,LOG
from . import __version__
from .translated import compact_models
from .identity import judge,background_evidence,profile_signal,query_coverage
from .recall import export_candidates
from .support import review_or_catalogue as review,audit_or_loci as audit
from .support.finalize_catalogue_review import export_unique


def repair_rank(m):
    intact=review.coding_intact(m) and m['canonical'] and m['has_start'] and m['has_stop']
    return (bool(m.get('qualified_reference_exact_matches')),bool(m.get('full_ref_support')) and intact,
        m.get('family')=='OR',profile_signal(m.get('or_profile_domains',[])),
        max((min(query_coverage(h),(h['se']-h['ss']+1)/h['slen']) for h in m.get('homology_hits',[])),default=0),
        m.get('or_bits',0),intact,-(max(e for s,e in m['exons'])-min(s for s,e in m['exons'])))


def eligible(m):
    # Weak homologs may enter REPAIR without being accepted as OR identities.
    from .tuning import policy
    t=policy(m)
    return (m['length']>=t['candidate_min_aa'] and (m.get('or_evalue',1)<=t['candidate_profile_evalue']
        or any(h['bits']>=t['candidate_bits'] and h['evalue']<=t['candidate_evalue'] and query_coverage(h)>=t['candidate_query_coverage'] for h in m.get('homology_hits',[]))
        or (m.get('alignment_reference_qualified',False) and m.get('ref_cov',0)>=.5 and m.get('identity',0)>=.15)
        or any(h.get('curated') and h['bits']>=30 and h['evalue']<=1e-3
               and query_coverage(h)>=.5 for h in m.get('homology_hits',[]))))


def model_key(m):
    return digest([m['contig'],m['strand'],m['exons'],m['cds']])


def plan_groups(groups,args,refs):
    from .engine import choose_refs
    runnable=[];unavailable=[]
    for m in groups:
        span=max(e for s,e in m['exons'])-min(s for s,e in m['exons'])+1
        reason='window_too_large' if span>args.max_window_bp else 'no_qualified_reference' if not choose_refs(m,refs,args.structure_references,args.reference_metadata,policy=getattr(args,'reconstruction_reference_policy','audited')) else ''
        if reason:unavailable.append(dict(model=m,reason=reason))
        else:runnable.append(m)
    # Scheduling only, not locus merging. Fixed anchors prevent a long bridge
    # from transitively joining tandem loci; every deferred model is retained.
    buckets=[]
    for m in runnable:
        match=next((b for b in buckets if b[0]['contig']==m['contig'] and b[0]['strand']==m['strand']
                    and max(review.overlap_fraction(b[0]['exons'],m['exons']),
                            review.overlap_fraction(m['exons'],b[0]['exons']))>=.8),None)
        if match is None:buckets.append([m])
        else:match.append(m)
    # Stronger alternate evidence contributes to the scheduling anchor, but
    # does not fuse coding paths or make an identity verdict. Short fragments
    # at already-visited loci cannot consume the whole budget before new loci.
    buckets.sort(key=lambda b:max(repair_rank(m) for m in b),reverse=True)
    scheduled=[b[i] for i in range(max(map(len,buckets),default=0)) for b in buckets if i<len(b)]
    return scheduled[:args.structure_max_loci],scheduled[args.structure_max_loci:],unavailable


def annotate(args,raw,refs,evidence,out):
    from .engine import structure,row_for,phylogeny,validate
    begin=time.monotonic();out=Path(out);folder=out/'annotation_qc';folder.mkdir(exist_ok=True)
    # raw is already expanded by orchestration. Nested alternatives are
    # reconstructed as flat groups below, never recursively copied per model.
    models=[copy.deepcopy({k:v for k,v in m.items() if k!='discovery_alternatives'}) for m in raw]
    params=SimpleNamespace(**vars(args));params.identity_policy='multichannel';params.force_structural_review=True
    params.flank_bp=args.structure_flank_bp
    old_policy=getattr(evidence.args,'identity_policy','legacy');evidence.args.identity_policy='multichannel'
    try:
        evidence.apply(models,'annotation_input_'+args.sample)
        groups=compact_models([m for m in models if eligible(m)],ranking=repair_rank)
        def assign(m):
            cid='ORA2_'+args.sample+'_'+digest([m['contig'],m['strand'],m['exons'],m['cds']])[:16]
            m.update(cid=cid,id=cid+'.original',source='original');return m
        groups=list(map(assign,groups));selected,deferred,unavailable=plan_groups(groups,args,refs)
        selected_ids={g['cid'] for g in selected};unavailable_ids={x['model']['cid'] for x in unavailable};group_keys={}
        for g in groups:
            disposition='reviewed_group' if g['cid'] in selected_ids else 'unavailable_window_or_reference' if g['cid'] in unavailable_ids else 'deferred_budget'
            for member in [g,*g.get('discovery_alternatives',[])]:group_keys[model_key(member)]=(g['cid'],disposition)
        trace={}
        for m in models:
            key=model_key(m);group,disposition=group_keys.get(key,('','low_priority_not_structurally_reviewed'))
            trace.setdefault(key,dict(input_candidate_id='ORC_'+digest((m['contig'],m['strand'],tuple(map(tuple,m['exons'])),m['cds']))[:20],
                repair_group_id=group,disposition=disposition,preliminary_identity=judge(m)['status']))
        review.write_table(folder/'input_disposition.tsv',list(trace.values()),['input_candidate_id','repair_group_id','disposition','preliminary_identity'])
        unreviewed=[m for m in models if trace[model_key(m)]['disposition']!='reviewed_group']
        export_candidates(unreviewed,out/'unreviewed_candidates')
        save(folder/'structure_plan.json',dict(input_models=len(models),eligible_models=sum(eligible(m) for m in models),groups=len(groups),
            selected=len(selected),deferred=len(deferred),unavailable=len(unavailable),force_review_complete_models=True,target_annotation_used=False))
        save(folder/'selected_models.json',selected);save(folder/'deferred_models.json',deferred)
        save(folder/'unavailable_groups.json',unavailable)
        LOG.info('Unified annotation: eligible groups=%d selected=%d deferred=%d',len(groups),len(selected),len(deferred))
        rows=[row_for(m) for m in selected];byid={};decisions={};tree_flags={};tree_status='not_needed'
        if rows:
            accepted,byid,decisions=structure(params,selected,out,evidence,refs)
            for row in rows:
                cid=row['candidate_id']
                if cid not in decisions:
                    original=next(m for m in byid[cid] if m['source']=='original')
                    damage=audit.damage_summary(cid,byid[cid],original,out/'structure/crosscheck')
                    decisions[cid]=audit.decide_locus(row,byid[cid],damage)
            LOG.info('[%s] Building the OR/GR diagnostic tree for %d model(s)', args.sample, len(rows))
            try:tree_flags=phylogeny(params,rows,decisions,refs,out);tree_status='completed'
            except Exception as error:
                LOG.exception('Diagnostic tree unavailable; no tree-only identity decisions')
                tree_status='failed';save(out/'phylogeny/diagnostic_failure.json',dict(error=str(error),identity_conservative=True))
        chosen=[decisions[r['candidate_id']]['model'] for r in rows]
        background_evidence(chosen,args,folder/'background')
        identity_rows=[];or_rows=[];or_decisions={};uncertain=[];excluded=[]
        for row in rows:
            cid=row['candidate_id'];d=decisions[cid];m=d['model'];flag=tree_flags.get(cid)
            m['identity_tree_evidence']=flag or {}
            m['identity_tree_reference_fingerprint']=evidence.config
            identity=judge(m,flag);status=identity['status']
            if status=='OR_supported':
                m['family']='OR';d['reason']+=';'+identity['reason']
                if d['category']=='complete' and (tree_status!='completed' or (flag or {}).get('identity_tree_vote')!='OR'):
                    d['confidence']='provisional_complete';d['reason']+=';identity_tree_not_reliable_for_reviewed_label'
                or_rows.append(row);or_decisions[cid]=d
            elif status.startswith('non_OR'):excluded.append(m)
            else:uncertain.append(m)
            identity_rows.append(dict(candidate_id=cid,selected_model_id=m['id'],identity=status,structure_category=d['category'],
                structure_confidence=d['confidence'],identity_reason=identity['reason'],structure_reason=d['reason'],**identity['signals']))
        review.write_table(folder/'identity_audit.tsv',identity_rows,list(identity_rows[0]) if identity_rows else ['candidate_id','identity'])
        # Preserve rejected/unresolved models; exclusion means exclusion from
        # the accepted OR catalogue, never erasure from the discovery data.
        export_candidates(uncertain,out/'identity_unresolved');export_candidates(excluded,out/'excluded_non_OR')
        if deferred:export_candidates(deferred,out/'structure_deferred')
        stats=dict(program='ORsearcher2',version=__version__,identity_counts=dict(Counter(r['identity'] for r in identity_rows)),tree_status=tree_status)
        if or_rows:
            review.export_results(or_rows,byid,or_decisions,SimpleNamespace(outdir=out),stats,{})
            unique=export_unique(or_rows,or_decisions,out,allow_containment=True)
            validation=validate(args.genome,out/'nonredundant')
        else:
            from .fasta import write_fasta
            final=out/'nonredundant';final.mkdir(exist_ok=True)
            review.write_table(final/'OR_catalogue.tsv',[],['candidate_id','sequence_id','strand','category','confidence'])
            (final/'OR_models.gff3').write_text('##gff-version 3\n')
            for kind in ['complete','partial','pseudogene']:
                write_fasta(final/f'OR_{kind}.proteins.faa',[]);write_fasta(final/f'OR_{kind}.cds.fna',[])
            unique=dict(representative_loci=0,category_counts={});validation=validate(args.genome,final)
        warnings=[]
        if deferred:warnings.append('structure_budget_deferred')
        if unavailable:warnings.append('structure_groups_without_window_or_reference')
        if tree_status=='failed':warnings.append('identity_tree_failed')
        rescuefile=out/'structure/rescue/summary.json'
        rescue=json.loads(rescuefile.read_text()) if rescuefile.exists() else dict(status='not_run')
        if rescue.get('deferred',0):warnings.append('bounded_structure_rescue_deferred')
        if rescue.get('status_counts',{}).get('timeout',0) or rescue.get('status_counts',{}).get('error',0):warnings.append('bounded_structure_rescue_unresolved_or_failed')
        if not getattr(args,'background_reference',None):warnings.append('non_OR_background_unavailable')
        failures=out/'structure/local_failures.tsv'
        if failures.exists() and len(failures.read_text().splitlines())>1:warnings.append('local_structure_failures')
        crossfile=out/'structure/crosscheck_summary.json'
        if crossfile.exists() and any(s['status']!='completed' for s in json.loads(crossfile.read_text())['statuses']):warnings.append('secondary_alignment_unresolved_or_failed')
        result=dict(status='completed_with_warnings' if warnings else 'completed',warnings=warnings,seconds=round(time.monotonic()-begin,2),
            selected_groups=len(selected),deferred_groups=len(deferred),unavailable_groups=len(unavailable),unreviewed_candidate_records=len({model_key(m) for m in unreviewed}),
            alternative_reference_rescue={k:v for k,v in rescue.items() if k!='statuses'},
            identity_counts=stats['identity_counts'],representatives=unique,validation=validation,
            phylogeny_status=tree_status,policy='Broad candidates retained; structure first, multichannel identity second; conflicting identities not deleted',
            limitations=['Computational support on one assembly, not experimental confirmation','Short/ambiguous models may remain unresolved','Genome-supported disruption can reflect assembly error; read/RNA validation not performed'])
        # Practical consolidation/repair is a final, separate policy. Preserve
        # the stricter evidence catalogue and expose only one primary results/.
        LOG.info('[%s] Final model selection, repair and export', args.sample)
        from .polish import polish
        polishing=polish(args,out,refs,evidence,groups=byid,decisions=or_decisions)
        result['polishing']=polishing
        if polishing['status']=='completed_with_warnings':
            result['status']='completed_with_warnings'
        save(folder/'summary.json',result)
        LOG.debug('Unified annotation details %s',result)
        LOG.info('[%s] Annotation finished | %s | results=%s',args.sample,result['status'],out/'results')
        return result
    finally:evidence.args.identity_policy=old_policy
