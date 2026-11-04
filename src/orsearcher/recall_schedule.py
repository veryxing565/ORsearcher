"""Annotation-blind local scheduling; never deletes discovery model records."""
from .tuning import policy as tuning_policy, reconstruction_metadata
from collections import defaultdict,Counter,deque


def origin(a):
    sources=set(a.get('origins',[a.get('origin','translated_hsp')]))
    return 'translated_hsp' if 'translated_hsp' in sources else 'metaeuk'


def rank(a):
    """Evidence priority, then source-specific score; not OR confirmation."""
    return (-bool(a.get('or_signal')), -a.get('or_bits',0),
            -a.get('chain_coverage',0),-a.get('distinct_references',0),
            -a.get('bits',0),a['contig'],a['strand'],a['start'],a['end'])


def exon_blocks(exons,max_span):
    """Separate long sparse sets at actual exon boundaries, not genome slicing.

    These are local search hints, NOT a decision to split a biological gene.
    Oversized individual exons remain explicit unschedulable hints.
    """
    groups=[]
    for s,e in sorted(exons):
        if not groups or e-groups[-1][0][0]+1>max_span:groups.append([(s,e)])
        else:groups[-1].append((s,e))
    return groups


def region(a,length,flank,max_window):
    for s,e in [(a.get('chain_start',a['start']),a.get('chain_end',a['end'])),(a['start'],a['end'])]:
        left=max(1,s-flank);right=min(length,e+flank)
        if left<=right and right-left+1<=max_window:return left,right
    return None


def reference_names(a,refs,metadata,limit):
    chosen=[];seen=set()
    for rid in a['references']:
        if rid not in refs or not metadata.get(rid,{}).get('reconstruction_allowed'):continue
        seq=refs[rid].sequence
        if seq in seen:continue
        chosen.append(rid);seen.add(seq)
        if len(chosen)>=limit:break
    return chosen


def spread_order(values):
    """Round-robin over 10kb bins; no hard cap or locus merging within a bin."""
    bins=defaultdict(list)
    for a in sorted(values,key=rank):bins[(a['contig'],a['strand'],a['start']//10000)].append(a)
    output=[]
    for depth in range(max((len(v) for v in bins.values()),default=0)):
        output.extend(sorted((v[depth] for v in bins.values() if len(v)>depth),key=rank))
    return output


def plan(args,hints,refs,lengths):
    """Budget counts runnable regions only. No GFF/annotation is consulted."""
    runnable=[];skipped=[]
    for a in hints:
        a=dict(a,references=list(a['references']))
        bounds=region(a,lengths[a['contig']],args.recall_flank_bp,args.recall_window_bp)
        names=reference_names(a,refs,reconstruction_metadata(args),args.recall_references)
        if bounds is None:skipped.append(dict(anchor=a,status='window_too_large'));continue
        if not names:skipped.append(dict(anchor=a,status='no_reference'));continue
        a.update(planned_start=bounds[0],planned_end=bounds[1],selected_references=names)
        runnable.append(a)
    # Equal opportunity between independent HSP and exon-set channels. Native
    # scores rank within each channel; they are never directly mixed globally.
    queues={k:deque(spread_order([a for a in runnable if origin(a)==k])) for k in ['translated_hsp','metaeuk']}
    ordered=[]
    while any(queues.values()):
        for k in ['translated_hsp','metaeuk']:
            if queues[k]:ordered.append(queues[k].popleft())
    selected=ordered[:args.recall_max_regions];deferred=ordered[len(selected):]
    summary=dict(input_hints=len(hints),runnable=len(runnable),selected=len(selected),deferred=len(deferred),
        skipped_counts=dict(Counter(x['status'] for x in skipped)),selected_origins=dict(Counter(origin(a) for a in selected)),
        policy='Runnable-only budget; source round-robin; spatial rounds without hard per-locus caps')
    return selected,deferred,skipped,summary
