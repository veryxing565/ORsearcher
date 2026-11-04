# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Whole-catalogue completeness, tandem-fusion and pseudogene audit.

Uses existing local Miniprot output INCLUDING nearby alignments, the full
second-aligner pass, OR/GR competition and the original genome. Predicted
annotation genes are contextual warnings, never accepted as ground truth.
"""
import argparse
import hashlib
import json
import logging
from collections import Counter,defaultdict
from pathlib import Path
from types import SimpleNamespace

from . import review_or_catalogue as review
from .finalize_catalogue_review import export_unique
from orsearcher.fasta import FastaRecord,read_fasta,write_fasta
from orsearcher.evidence import ga_thresholds, parse_domains, profile_cutoff
from orsearcher.models import parse_miniprot_gff,attach_reference_lengths,parse_attributes
from orsearcher.sharding import Segment

LOG=logging.getLogger('locus_audit')


def bases(exons):return sum(e-s+1 for s,e in exons)


def shared(a,b):
    return sum(max(0,min(e,v)-max(s,u)+1) for s,e in a for u,v in b)


def bounds(model):return min(s for s,e in model['exons']),max(e for s,e in model['exons'])


def tool(model):
    # Global original and local Miniprot predictions are the SAME aligner.
    # Carry the true aligner when discovery used Exonerate and orchestration
    # placed that model into its generic 'original' slot. It is not Miniprot.
    return model.get('alignment_tool') or ('exonerate' if 'exonerate' in model['source'] else 'miniprot')


def same_region(model,original):
    a,b=bounds(model);c,d=bounds(original)
    overlap=max(0,min(b,d)-max(a,c)+1)/min(b-a+1,d-c+1)
    # Existing CDS may be wrong. Conservatively permit a replacement with
    # substantial shared coding AND approximately the same genomic boundaries.
    return model['overlap']>=.8 or model['overlap']>=.5 and overlap>=.8 and abs(a-c)<=300 and abs(b-d)<=300


def near_target(model,original):
    a,b=bounds(model);c,d=bounds(original)
    region_overlap=max(0,min(b,d)-max(a,c)+1)/min(b-a+1,d-c+1)
    return model['overlap']>=.5 and region_overlap>=.8


def junctions(model):
    exons=sorted(model['exons'])
    return {(left[1]+1,right[0]-1) for left,right in zip(exons,exons[1:])}


def structurally_agrees(a,b):
    if min(review.overlap_fraction(a['exons'],b['exons']),review.overlap_fraction(b['exons'],a['exons']))<.95:return False
    ja,jb=junctions(a),junctions(b)
    if not ja and not jb:return True
    return ja==jb


def endpoints_agree(a,b):
    return all(abs(x-y)<=3 for x,y in zip(bounds(a),bounds(b)))

def external_alignment(model):
    ref=model.get('ref','')
    return bool(ref and not ref.startswith(('adaptive_','annotation_')) and model.get('alignment_reference_qualified',True))

def external_identity(model):
    return model.get('identity',0) if external_alignment(model) else 0.

def reference_key(model):
    return model.get('reference_sequence_key') or model.get('ref','')

def external_exact(model):
    return bool(model.get('qualified_reference_exact_matches',model.get('reference_exact_matches',[])))

def compatible_support(chosen,models,independent=False):
    """One vote per real tool/reference SEQUENCE, not per terminal variant."""
    exact=set(chosen.get('reference_exact_matches',[])) if independent else set()
    unique={}
    for m in models:
        if not external_alignment(m) or m['ref'] in exact or m.get('ref_cov',0)<.8:continue
        if not(review.coding_intact(m) and m['canonical'] and m['family']=='OR'
               and structurally_agrees(m,chosen) and endpoints_agree(m,chosen)):continue
        unique[(tool(m),reference_key(m))]=m
    return list(unique.values())

def selection_rank(model,models):
    """External compatibility first; no length or self-identity tie breaker."""
    support=compatible_support(model,models)
    other_refs={reference_key(m) for m in support if tool(m)!=tool(model)} if external_alignment(model) else set()
    endpoint_error=sum(abs(model.get('terminal',{}).get(k,0)) for k in ['n_difference','c_difference'])
    # A larger homolog list can reflect uneven database sampling or a shorter
    # query, not better endpoints. Count independent STRUCTURAL references
    # above; whole-protein support is presence/absence rather than a vote.
    return (external_exact(model),external_alignment(model) and model.get('ref_cov',0)>=.8,
            bool(model.get('full_ref_support',[])),len(other_refs),
            -endpoint_error,external_identity(model)*min(model.get('ref_cov',0),1),model.get('or_bits',0))


def structure_candidate(model):
    return model.get('family')=='OR' or model.get('structure_or_eligible',False)


def structural_homology_hits(model):
    """Comparative evidence, NOT identity confirmation or repair admission.

    Query fractions must not erase the same local HSP on a longer coding
    path. Input homology_hits already underwent reference/domain audit.
    """
    return [h for h in model.get('homology_hits',[])
            if h['bits']>=50 and h['evalue']<=1e-5 and h['identity']>=15
            and h['se']-h['ss']+1>=80 and h['slen']>0]


def partial_selection_evidence(model,models):
    hits=[h for m in models for h in structural_homology_hits(m)]
    key=lambda h:h.get('reference_sequence_key',h['ref'])
    anchor=max(hits,key=lambda h:(h['bits'],h.get('curated',False),key(h)),default=None)
    shared=[h for h in structural_homology_hits(model) if anchor is not None and key(h)==key(anchor)]
    coverage=max(((h['se']-h['ss']+1)/max(1,h['slen']) for h in shared),default=0)
    bits=max((h['bits'] for h in shared),default=0)
    strongest=max(shared,key=lambda h:(h['bits'],h['identity']),default=None)
    return dict(comparison_reference=anchor['ref'] if anchor else '',
        comparison_reference_key=key(anchor) if anchor else '',
        reference_coverage=min(coverage,1),shared_bits=bits,
        shared_identity=strongest['identity'] if strongest else 0,
        query_fraction_diagnostic=(strongest['qe']-strongest['qs']+1)/max(1,strongest['qlen']) if strongest else 0,
        query_fraction_used_in_rank=False,
        coding_endpoints=int(bool(model['has_start']))+int(bool(model['has_stop'])),
        original_anchor=model['source']=='original')


def partial_selection_rank(model,models):
    """Same-reference support; no query denominator or mixed native identity.

    Equal evidence retains coding endpoints/valid original provenance, not the
    longest peptide. Damaged/noncanonical/fused alternatives remain gated.
    """
    support=partial_selection_evidence(model,models)
    endpoint_error=sum(abs(model.get('terminal',{}).get(k,0)) for k in ['n_difference','c_difference'])
    return (review.coding_intact(model),external_exact(model),int((support['reference_coverage']+1e-9)*20),
            support['shared_bits'],support['shared_identity'],support['coding_endpoints'],
            support['original_anchor'],-endpoint_error,model.get('or_bits',0),model['id'])

def internal_splice_conflicts(chosen,models):
    """Comparable external splice paths stay unresolved, not falsely complete."""
    if external_exact(chosen):return []
    support=compatible_support(chosen,models,independent=True)
    if {'miniprot','exonerate'}<={tool(m) for m in support}:return []
    return [m for m in models if m is not chosen and external_alignment(m) and m['ref_cov']>=.8
            and review.coding_intact(m) and m['canonical'] and m['family']=='OR' and not m['fusion_pairs']
            and min(review.overlap_fraction(chosen['exons'],m['exons']),review.overlap_fraction(m['exons'],chosen['exons']))>=.8
            and junctions(m)!=junctions(chosen) and endpoints_agree(m,chosen)
            and m.get('or_bits',0)>=chosen.get('or_bits',0)*.95
            and abs(external_identity(m)-external_identity(chosen))<=.02]


def truncation_conflicts(chosen,models):
    """Keep a plausible longer same-locus path as unresolved, not silently lost."""
    return [m for m in models if m is not chosen and near_target(m,chosen)
            and m['family']=='OR' and m['canonical'] and not m['fusion_pairs']
            and review.coding_intact(m) and m.get('or_ga',False)
            and external_alignment(m) and m['ref_cov']>=.8 and external_identity(m)>=max(.3,external_identity(chosen)-.1)
            and m['length']>=chosen['length']+max(40,int(chosen['length']*.1))
            and review.overlap_fraction(chosen['exons'],m['exons'])>=.9]


def collect_neighbors(rows,run,refs):
    unique={}
    for row in rows:
        cid=row['candidate_id'];case=run/'cases'/cid
        if not all((case/p).is_file() for p in ['region.fna','region.json','miniprot.gff3']):continue
        seq=next(read_fasta(case/'region.fna')).sequence
        meta=json.loads((case/'region.json').read_text());offset=meta['offset']
        segment=Segment('region',row['sequence_id'],'',offset,offset+len(seq),offset+len(seq),'region.fna',False)
        parsed=parse_miniprot_gff(case/'miniprot.gff3',{'region':segment})
        attach_reference_lengths(parsed,{r:len(v.sequence) for r,v in refs.items()})
        for item in parsed:
            features=sorted([f for f in item.features if f.feature_type=='CDS'],key=lambda f:f.start,reverse=item.strand=='-')
            if not features or item.coverage<.7:continue
            try:m=review.build_model(cid,'neighbor_miniprot',seq,offset,[(f.start,f.end) for f in features],item.strand,
                  int(features[0].phase) if features[0].phase.isdigit() else 0,item.frameshifts,item.reference_id,
                  item.reference_start,item.reference_end,item.reference_length,item.identity)
            except ValueError:continue
            m['alignment']=dict(qstart=item.reference_start,qend=item.reference_end,qlen=item.reference_length)
            m['contig']=row['sequence_id']
            variants=[m,*review.terminal_variants(m,seq,offset,refs)]
            for v in variants:
                v['contig']=row['sequence_id']
                key=json.dumps([v['contig'],v['strand'],v['exons'],v['ref']],sort_keys=True)
                v['id']='NEIGH_'+hashlib.sha256(key.encode()).hexdigest()[:18]
                unique[v['id']]=v
    return list(unique.values())


def anchor_groups(models):
    good=[m for m in models if m['family']=='OR' and m['canonical'] and not m['frameshifts']
          and len(m['stops'])<=1 and not m['ambiguous'] and not m['remainder']
          and .8<=m['ref_cov']<=1.05 and m['identity']>=.3 and 250<=m['length']<=650
          and (m['or_ga'] or m['or_bits']>=150)]
    good.sort(key=lambda m:(m['identity'],m['ref_cov'],m['or_bits'],-(bounds(m)[1]-bounds(m)[0])),reverse=True)
    groups=[]
    for m in good:
        match=None
        for g in groups:
            anchor=g['model']
            if anchor['contig']==m['contig'] and anchor['strand']==m['strand'] and min(
               review.overlap_fraction(anchor['exons'],m['exons']),review.overlap_fraction(m['exons'],anchor['exons']))>=.65:
                match=g;break
        if match is None:
            groups.append(dict(id='ORLOC_'+str(len(groups)+1).zfill(4),model=m,members=[m]))
        else:match['members'].append(m)
    return groups


def fusion_links(model,anchors):
    touches=[]
    for g in anchors:
        a=g['model']
        if a['contig']!=model['contig'] or a['strand']!=model['strand']:continue
        overlap=shared(model['exons'],a['exons'])
        if overlap>=120 and overlap/bases(model['exons'])>=.15 and overlap/bases(a['exons'])>=.15:
            touches.append(g)
    pairs=[]
    for i,left in enumerate(touches):
        for right in touches[i+1:]:
            ls,le=bounds(left['model']);rs,re=bounds(right['model'])
            # A long intron spanning another gene is not by itself a fusion.
            # Actual candidate CDS must use BOTH independent spatial loci.
            if le<rs or re<ls:pairs.append((left['id'],right['id']))
    return pairs


def protein_repeats(models,args):
    evidence=args.outdir/'protein_structure';evidence.mkdir(exist_ok=True)
    fasta=evidence/'models.faa'
    write_fasta(fasta,[FastaRecord(m['id'],'',m['protein'].replace('*','X')) for m in models if m['protein']])
    review.command([args.hmmsearch,'--cpu',str(args.workers),'--noali','--domtblout',str(evidence/'or.domtbl'),str(args.or_hmm),str(fasta)],evidence/'or.out',180)
    domains=defaultdict(list)
    cutoffs=ga_thresholds(args.or_hmm)
    for mid, hits in parse_domains(evidence/'or.domtbl').items():
        for hit in hits:
            cutoff=profile_cutoff(hit,cutoffs)
            if hit['domain_evalue']<=1e-3 and hit['coverage']>=.5 and (cutoff is None or hit['domain_score']>=cutoff[1]):
                domains[mid].append((hit['start'],hit['end'],hit['profile']))
    # Only selected/alternative TARGET models use multiple HSP mode. Nearby
    # models are useful spatial evidence but do not need expensive repeat scans.
    target=[m for m in models if not m['source'].startswith('neighbor')]
    write_fasta(evidence/'target_models.faa',[FastaRecord(m['id'],'',m['protein'].replace('*','X')) for m in target if m['protein']])
    path=evidence/'or_multi_hsp.tsv'
    review.command([args.diamond,'blastp','--query',str(evidence/'target_models.faa'),'--db',str(args.or_db),'--out',str(path),
        '--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore',
        '--max-hsps','3','--max-target-seqs','5','--sensitive','--evalue','1e-5','--threads',str(args.workers),'--quiet'],timeout=180)
    hits=defaultdict(list)
    for line in path.read_text().splitlines():
        f=line.split('\t')
        if int(f[3])>=80 and float(f[2])>=25 and float(f[11])>=40:
            hits[(f[0],f[1])].append((int(f[6]),int(f[7]),int(f[8]),int(f[9])))
    duplicated=set()
    for (mid,rid),hsps in hits.items():
        for i,a in enumerate(hsps):
            for b in hsps[i+1:]:
                qo=max(0,min(a[1],b[1])-max(a[0],b[0])+1)
                so=max(0,min(a[3],b[3])-max(a[2],b[2])+1)
                if qo<=5 and so/min(a[3]-a[2]+1,b[3]-b[2]+1)>=.6:duplicated.add(mid)
    for m in models:
        ds=domains[m['id']]
        m['repeated_or_domains']=any(a[2]==b[2] and min(a[1],b[1])-max(a[0],b[0])+1<=5
                                     for i,a in enumerate(ds) for b in ds[i+1:])
        m['repeated_reference_regions']=m['id'] in duplicated


def annotation_context(path,models):
    genes={}
    for line in path.read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split('\t')
        if len(f)!=9 or f[2]!='CDS':continue
        parent=parse_attributes(f[8]).get('Parent','')
        key=(f[0],f[6],parent)
        genes.setdefault(key,[]).append((int(f[3]),int(f[4])))
    bycontig=defaultdict(list)
    for (seqid,strand,gid),exons in genes.items():bycontig[(seqid,strand)].append((gid,exons))
    for m in models:
        ms,me=bounds(m);overlap=[]
        for gid,exons in bycontig[(m['contig'],m['strand'])]:
            if max(e for s,e in exons)<ms or min(s for s,e in exons)>me:continue
            bp=shared(m['exons'],exons)
            if bp>=120 and bp/bases(m['exons'])>=.15:overlap.append(gid)
        m['annotation_gene_ids']=overlap


def damage_summary(cid,models,original,crosscheck,frame_cache=None):
    trusted=[m for m in models if near_target(m,original) and m['family']=='OR' and m['canonical']
             and external_alignment(m) and m['ref_cov']>=.75 and m['identity']>=.25 and not m['ambiguous'] and not m['fusion_pairs']]
    stops=defaultdict(list)
    for m in trusted:
        if m['frameshifts']:continue
        for s in m['stops']:
            if s['edge_distance']>=30 and s['ambiguous_context']==0:
                stops[tuple(s['coordinates'])].append(m)
    strongest=None
    for coords,support in stops.items():
        refs={m['ref'] for m in support};tools={tool(m) for m in support}
        independent={m.get('reference_sequence_key') or m['ref'] for m in support}
        score=(len(independent),len(tools))
        if strongest is None or score>strongest['score']:
            strongest=dict(kind='premature_stop',coordinates=list(coords),references=sorted(refs),independent_reference_count=len(independent),tools=sorted(tools),score=score,models=support)
    events=[];case=crosscheck/'cases'/cid
    batches=frame_cache if frame_cache is not None else [
          (aligner,json.loads((case/filename).read_text())) for filename,aligner in
          [('miniprot_frame_events.json','miniprot'),('exonerate_frame_events.json','exonerate')] if (case/filename).exists()]
    for aligner,alignments in batches:
        for a in alignments:
            if a['strand']!=original['strand']:continue
            matching=[m for m in trusted if tool(m)==aligner and m['ref']==a['reference']
                      and abs(bounds(m)[0]-a['start'])<15 and abs(bounds(m)[1]-a['end'])<15 and m['ref_cov']>=.8 and m['identity']>=.5]
            if not matching:continue
            for e in a['events']:
                coords=e['coordinates'] or [e.get('boundary',0)]
                inside=[m for m in matching if any(s+15<=min(coords)<=max(coords)<=v-15 for s,v in m['exons'])]
                if inside and e['context_ambiguous']==0:
                    events.append(dict(position=min(coords),coordinates=coords,reference=a['reference'],tool=aligner,models=inside))
    groups=[]
    for e in sorted(events,key=lambda x:x['position']):
        g=next((g for g in groups if abs(g[0]['position']-e['position'])<=2),None)
        if g is None:groups.append([e])
        else:g.append(e)
    frame=None
    for g in groups:
        refs={e['reference'] for e in g};tools={e['tool'] for e in g}
        independent={m.get('reference_sequence_key') or m['ref'] for e in g for m in e['models']}
        score=(len(tools),len(independent))
        if frame is None or score>frame['score']:
            frame=dict(kind='frame_disruption',coordinates=g[0]['coordinates'],references=sorted(refs),independent_reference_count=len(independent),tools=sorted(tools),score=score,
                       models=[m for e in g for m in e['models']])
    if strongest and strongest['independent_reference_count']>=3 and {'miniprot','exonerate'}<=set(strongest['tools']):return strongest
    if frame and frame['independent_reference_count']>=2 and {'miniprot','exonerate'}<=set(frame['tools']):return frame
    return strongest or frame


def decide_locus(row,models,damage):
    original=next(m for m in models if m['source']=='original')
    eligible=[m for m in models if same_region(m,original) and not m['fusion_pairs'] and not (m['repeated_or_domains'] and m['repeated_reference_regions'])]
    good=[]
    for m in eligible:
        # Use the normal joint completeness criteria, but a genome-backed same-
        # locus reconstruction may replace incorrect internal exons of the old model.
        proxy=dict(m,overlap=1.,bridge=False)
        if review.complete_model(proxy):good.append(m)
    credible_intact=[m for m in models if near_target(m,original) and review.coding_intact(m) and m['canonical'] and structure_candidate(m)
                    and m['ref_cov']>=.8 and m['identity']>=.3 and m['length']>=250 and not m['fusion_pairs']]
    strong_damage=damage and {'miniprot','exonerate'}<=set(damage['tools']) and damage.get('independent_reference_count',len(damage['references']))>=(3 if damage['kind']=='premature_stop' else 2)
    # A distant-reference intact splice path must not silently override a much
    # closer full-coverage disrupted path. Even a frameshift near a splice edge
    # (not sufficient to CONFIRM a pseudogene) makes that structure unresolved.
    closest_disrupted=[m for m in models if near_target(m,original) and m['family']=='OR' and m['canonical']
                       and external_alignment(m) and m['ref_cov']>=.9 and m['identity']>=.85 and not m['ambiguous']
                       and (m['frameshifts'] or m['stops']) and not m['fusion_pairs']]
    if good and closest_disrupted and max(m['identity'] for m in closest_disrupted)>=max(external_identity(m) for m in good)+.2:
        chosen=max(good,key=lambda m:selection_rank(m,good))
        return dict(category='partial',confidence='unresolved',model=chosen,reason='closer_disrupted_homolog_conflicts_with_distant_intact_splice_model',
                    stop_reference_count=len(damage['references']) if damage else 0,
                    stop_coordinates=','.join(map(str,damage['coordinates'])) if damage else '',crosscheck_sources=','.join(damage['tools']) if damage else '')
    if good and strong_damage:
        # A plausible intact and a reproducibly disrupted splice model at the
        # same locus are a real structural conflict, NOT proof either is correct.
        chosen=max(good,key=lambda m:selection_rank(m,good))
        return dict(category='partial',confidence='unresolved',model=chosen,reason='intact_and_disrupted_same_locus_models_conflict',
                    stop_reference_count=len(damage['references']),stop_coordinates=','.join(map(str,damage['coordinates'])),crosscheck_sources=','.join(damage['tools']))
    if good:
        ordered=sorted(good,key=lambda m:(str(m['exons']),m.get('ref',''),tool(m),m['id']))
        chosen=max(ordered,key=lambda m:selection_rank(m,good))
        longer=truncation_conflicts(chosen,models)
        if longer:
            return dict(category='partial',confidence='unresolved',model=chosen,reason='longer_same_locus_OR_model_conflicts_with_complete_endpoints',
                        stop_reference_count=0,stop_coordinates='',crosscheck_sources='')
        conflicts=internal_splice_conflicts(chosen,good)
        if conflicts:
            return dict(category='partial',confidence='unresolved',model=chosen,reason='comparable_external_internal_splice_paths_conflict',
                        stop_reference_count=0,stop_coordinates='',crosscheck_sources='')
        supporting=compatible_support(chosen,eligible,independent=True)
        tools={tool(m) for m in supporting};refs={m['ref'] for m in supporting if m['ref']}
        confidence='reviewed_complete' if {'miniprot','exonerate'}<=tools and refs else 'provisional_complete'
        return dict(category='complete',confidence=confidence,model=chosen,reason='complete_coding_model_with_whole_catalogue_locus_checks',
                    stop_reference_count=0,stop_coordinates='',crosscheck_sources=','.join(sorted(tools)))
    if strong_damage and not credible_intact:
        chosen=max(damage['models'],key=lambda m:(m['ref_cov'],m['identity'],-m['frameshifts'],-len(m['stops'])))
        return dict(category='pseudogene',confidence='genome_supported',model=chosen,reason='reproducible_exonic_disruption_no_credible_intact_alternative',
             stop_reference_count=len(damage['references']),stop_coordinates=','.join(map(str,damage['coordinates'])),crosscheck_sources=','.join(damage['tools']))
    viable=[m for m in models if near_target(m,original) and structure_candidate(m) and not m['fusion_pairs'] and m['canonical']
            and not(m.get('repeated_or_domains') and m.get('repeated_reference_regions'))]
    chosen=max(viable,key=lambda m:partial_selection_rank(m,viable)) if viable else original
    return dict(category='partial',confidence='unresolved',model=chosen,reason='incomplete_or_conflicting_structure_or_insufficient_disruption_evidence',
                stop_reference_count=len(damage['references']) if damage else 0,
                stop_coordinates=','.join(map(str,damage['coordinates'])) if damage else '',crosscheck_sources=','.join(damage['tools']) if damage else '')


