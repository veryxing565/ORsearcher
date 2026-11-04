"""Conservative computational identity, separate from coding completeness.

No score/length/tree condition alone establishes or deletes an OR. Conflicts
and missing evidence remain explicit. Thresholds are development defaults,
not experimentally validated probabilities or independent biological tests.
"""
from pathlib import Path
from collections import defaultdict
from .runtime import execute,save,file_digest,digest
from .fasta import FastaRecord,write_fasta


def query_coverage(hit):
    return (hit['qe']-hit['qs']+1)/max(1,hit['qlen'])


def credible_hits(hits,t=None):
    from .tuning import DEFAULTS
    t=t or DEFAULTS
    return [h for h in hits if h['evalue']<=t['identity_evalue'] and h['bits']>=t['identity_bits']
            and query_coverage(h)>=t['identity_query_coverage'] and h['identity']>=t['identity_min_identity']]


def profile_signal(domains,t=None):
    from .tuning import DEFAULTS
    t=t or DEFAULTS
    return any(d['domain_evalue']<=t['profile_evalue'] and d['coverage']>=t['profile_coverage'] for d in domains)


def distant_hits(hits,t=None):
    """Not sufficient alone: used only alongside an OR profile/consensus."""
    from .tuning import DEFAULTS
    t=t or DEFAULTS
    return [h for h in hits if h['evalue']<=t['distant_evalue'] and h['bits']>=t['distant_bits']
            and query_coverage(h)>=t['distant_query_coverage'] and h['identity']>=t['distant_min_identity']]


def full_curated_consensus(hits):
    """Near-whole strong curated homolog PLUS a distinct significant homolog.

    Profile-negative ORs need not depend on one diagnostic tree. A lone strong
    hit, duplicated references or a local high-identity hit is insufficient.
    """
    distant=distant_hits(hits)
    whole=[h for h in distant if h.get('curated') and h['bits']>=150
           and h['evalue']<=1e-20 and h['identity']>=35 and query_coverage(h)>=.85
           and (h['se']-h['ss']+1)/max(1,h['slen'])>=.85]
    secondary={h.get('reference_sequence_key',h['ref']) for h in distant
               if h.get('curated') and h['bits']>=45 and h['evalue']<=1e-5}
    return any(secondary-{h.get('reference_sequence_key',h['ref'])} for h in whole)


def judge(model,tree=None):
    """Tree votes require reference separation AND local/occupancy QC."""
    from .tuning import policy
    t=policy(model)
    if tree is None and t['identity_use_tree'] and model.get('identity_tree_reference_fingerprint')==model.get('_evidence_config') and model.get('_evidence_config') is not None:tree=model.get('identity_tree_evidence')
    oh=credible_hits(model.get('homology_hits',[]),t);gh=credible_hits(model.get('gr_homology_hits',[]),t)
    od=distant_hits(model.get('homology_hits',[]),t);gd=distant_hits(model.get('gr_homology_hits',[]),t)
    # Competing family scores use the SAME permissive metrics. Discarding a
    # distant GR competitor while accepting a distant OR would be asymmetric.
    ob=max((h['bits'] for h in od),default=0);gb=max((h['bits'] for h in gd),default=0)
    op=profile_signal(model.get('or_profile_domains',[]),t)
    gp=any(profile_signal(ds,t) for ds in model.get('gr_profile_domains',[]))
    keys={h.get('reference_sequence_key',h['ref']) for h in oh if h.get('curated')}
    exact=bool(model.get('qualified_reference_exact_matches'))
    vote=(tree or {}).get('identity_tree_vote','unresolved')
    consensus=len(keys)>=2 and ob>=100
    positive=bool(oh) and ob-gb>=t['family_score_margin'] and (op or consensus or exact or vote=='OR')
    strong_profile=any(d['domain_evalue']<=1e-5 and d['coverage']>=.35
                       for d in model.get('or_profile_domains',[]))
    distant_keys={h.get('reference_sequence_key',h['ref']) for h in od if h.get('curated')}
    remote_full=any(h.get('curated') and h['bits']>=50 and h['evalue']<=1e-8
                    and query_coverage(h)>=.5 for h in od)
    remote_consensus=len(distant_keys)>=2 and any(h['bits']>=40 and h['evalue']<=1e-5 for h in od)
    distant_positive=strong_profile and ob-gb>=t['distant_score_margin'] and (remote_full or remote_consensus)
    profile_independent=full_curated_consensus(model.get('homology_hits',[])) and ob-gb>=30 and not gp
    positive=positive or distant_positive or profile_independent
    negative=bool(gh) and gb-ob>=t['family_score_margin'] and (gp or vote=='GR')
    background=[h for h in model.get('background_hits',[]) if h['evalue']<=1e-15 and h['bits']>=100
        and h['identity']>=70 and query_coverage(h)>=.85
        and (h['se']-h['ss']+1)/max(1,h['slen'])>=.5]
    bb=max((h['bits'] for h in background),default=0)
    other=bool(background) and bb-max(ob,gb)>=t['background_score_margin']
    signals=dict(or_profile=op,gr_profile=gp,or_homology_bits=ob,gr_homology_bits=gb,
        curated_OR_sequences=len(keys),qualified_exact_OR=exact,tree_vote=vote,background_bits=bb,
        distant_OR_sequences=len(distant_keys),distant_OR_path=distant_positive,
        profile_independent_full_consensus=profile_independent,
        tree_disagreement=positive and vote=='GR')
    # A tree alone never vetoes a multi-channel OR. Real GR/background
    # competition still causes a conflict; the tree disagreement is audited.
    if (positive and (negative or other or (vote=='GR' and gb>=ob-t['family_score_margin'] and gb>0))) or (negative and (op or vote=='OR')) or (other and op):
        return dict(status='identity_unresolved',reason='conflicting_OR_GR_background_or_tree_evidence',signals=signals)
    if positive:return dict(status='OR_supported',reason='OR_supported_tree_disagreement_without_GR_competition' if vote=='GR'
        else 'near_full_curated_OR_with_distinct_homolog_support' if profile_independent and not op
        else 'distant_OR_homology_with_strong_profile_and_curated_support' if distant_positive and not oh
        else 'OR_homology_with_profile_curated_consensus_exact_or_reliable_tree_support',signals=signals)
    if negative:return dict(status='non_OR_GR_supported',reason='GR_competition_with_profile_or_reliable_tree_no_OR_conflict',signals=signals)
    if other and not op and not positive:
        return dict(status='non_OR_other_supported',reason='near_full_curated_non_OR_match_with_large_competitive_margin_no_OR_profile_support',signals=signals)
    return dict(status='identity_unresolved',reason='insufficient_or_unbalanced_family_evidence',signals=signals)


def structural_family(model):
    result=judge(model)
    model['legacy_family']=model.get('family','unknown');model['preliminary_identity']=result
    model['family']='OR' if result['status']=='OR_supported' else 'GR_likely' if result['status']=='non_OR_GR_supported' else 'OR_GR_uncertain'
    # Structural eligibility is not a final identity verdict. Weak significant
    # OR homologs may have their coding structures checked before tree evidence.
    from .tuning import policy
    t=policy(model)
    or_hits=distant_hits(model.get('homology_hits',[]),t);gr_hits=distant_hits(model.get('gr_homology_hits',[]),t)
    ob=max((h['bits'] for h in or_hits),default=0);gb=max((h['bits'] for h in gr_hits),default=0)
    model['structure_or_eligible']=bool(or_hits) and ob-gb>=t['family_score_margin'] and not result['status'].startswith('non_OR')
    if result['signals']['profile_independent_full_consensus']:
        model['family_support']='curated_full_homology_consensus'


def background_evidence(models,args,folder):
    """A bounded protein-only check, never searches genomes against this DB."""
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    reference=getattr(args,'background_reference',None)
    if reference is None:
        save(folder/'summary.json',dict(status='disabled',limitation='Other-protein exclusion unavailable; unknown candidates retained'))
        return
    reference=Path(reference);db=folder/'negative';query=folder/'proteins.faa'
    sequences={digest(m['protein'].replace('*','X')):m['protein'].replace('*','X') for m in models if m['protein']}
    write_fasta(query,[FastaRecord(k,'',v) for k,v in sequences.items()])
    if not sequences:return
    execute(['diamond','makedb','--in',reference,'--db',db,'--threads',args.threads],folder/'makedb.log',inputs=[reference])
    output=folder/'hits.tsv'
    execute(['diamond','blastp','--query',query,'--db',db,'--out',output,'--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore',
        '--max-target-seqs','5','--sensitive','--evalue','1e-5','--threads',args.threads,'--quiet'],folder/'search.log',inputs=[query,Path(str(db)+'.dmnd')])
    hits=defaultdict(list)
    for line in output.read_text().splitlines():
        f=line.split('\t');hits[f[0]].append(dict(ref=f[1],identity=float(f[2]),length=int(f[3]),qlen=int(f[4]),slen=int(f[5]),
            qs=int(f[6]),qe=int(f[7]),ss=int(f[8]),se=int(f[9]),evalue=float(f[10]),bits=float(f[11])))
    for m in models:m['background_hits']=hits[digest(m['protein'].replace('*','X'))]
    save(folder/'summary.json',dict(status='completed',reference=str(reference),reference_sha256=file_digest(reference),unique_queries=len(sequences),
        policy='Protein-only reviewed non-OR competition; no deletion on a hit alone; strong OR conflicts retained'))
