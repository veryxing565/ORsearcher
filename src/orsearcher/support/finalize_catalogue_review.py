# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Reassess cached genomic models, and export conservative same-locus aliases.

Writes a separate output; raw results and all alternative models are preserved.
Never merges by protein identity alone. Anchor-based assignment avoids transitive
merging through a model that overlaps two distinct loci.
"""
import argparse
import json
import logging
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

from . import review_or_catalogue as review
from orsearcher.fasta import FastaRecord,read_fasta,write_fasta


def same_locus(left,right):
    a=review.overlap_fraction(left["exons"],right["exons"])
    b=review.overlap_fraction(right["exons"],left["exons"])
    if min(a,b)>=.8:return True
    if min(a,b)<.65:return False
    def junctions(m):
        exons=sorted(m['exons'])
        return {(x[1]+1,y[0]-1) for x,y in zip(exons,exons[1:])}
    ja,jb=junctions(left),junctions(right)
    # Shared coding and splice structure, not mere neighbouring proximity.
    if ja and jb:return len(ja&jb)/max(len(ja),len(jb))>=.6
    if not ja and not jb:return (left['exons'][0][0]-right['exons'][0][0])%3==0
    return False

def shared_reading_frame(left,right):
    """Same codon offset on every shared CDS block (both strands)."""
    if left['strand']!=right['strand']:return False
    def blocks(m):
        cumulative=0;out=[]
        for s,e in sorted(m['exons'],reverse=m['strand']=='-'):
            out.append((s,e,cumulative));cumulative+=e-s+1
        return out
    found=False
    for s,e,offset in blocks(left):
        for t,u,other in blocks(right):
            if max(s,t)>min(e,u):continue
            found=True;base=max(s,t) if left['strand']=='+' else min(e,u)
            a=offset+(base-s if left['strand']=='+' else e-base)
            b=other+(base-t if right['strand']=='+' else u-base)
            if (a-b)%3:return False
    return found

def terminal_equivalent(left,right):
    """Technical endpoint alternatives, NOT a protein-similarity locus merge."""
    if not all(review.coding_intact(m) and m['canonical'] and m.get('family')=='OR'
               and not m.get('fusion_pairs') and not m.get('repeated_or_domains') for m in [left,right]):return False
    if not shared_reading_frame(left,right):return False
    a=review.overlap_fraction(left['exons'],right['exons']);b=review.overlap_fraction(right['exons'],left['exons'])
    if min(a,b)<.98:return False
    le,re=sorted(left['exons']),sorted(right['exons'])
    if [(x[1],y[0]) for x,y in zip(le,le[1:])]!=[(x[1],y[0]) for x,y in zip(re,re[1:])]:return False
    return abs(le[0][0]-re[0][0])<=12 and abs(le[-1][1]-re[-1][1])<=12

def splice_compatible(left,right):
    """Neither fragment translates the other's skipped intron as coding DNA."""
    for a,b in [(left,right),(right,left)]:
        exons=sorted(a['exons'])
        for x,y in zip(exons,exons[1:]):
            if any(max(0,min(y[0]-1,e)-max(x[1]+1,s)+1) for s,e in b['exons']):return False
    return True

def fragment_alias(left,right):
    """Near-contained intact partials only; never connect disjoint fragments."""
    if not all(review.coding_intact(m) and m['canonical'] and m.get('family')=='OR'
               and not m.get('fusion_pairs') and not m.get('repeated_or_domains') for m in [left,right]):return False
    if not shared_reading_frame(left,right) or not splice_compatible(left,right):return False
    lengths=[sum(e-s+1 for s,e in m['exons']) for m in [left,right]]
    small,big=(left,right) if lengths[0]<=lengths[1] else (right,left)
    a=review.overlap_fraction(small['exons'],big['exons']);b=review.overlap_fraction(big['exons'],small['exons'])
    outside=(1-a)*min(lengths)
    return a>=.9 and b>=.55 and outside<=30.000001


def export_unique(rows,decisions,outdir,allow_containment=False):
    """Only reciprocally overlapping models of the same category are aliased.

    Identical-coordinate duplicate models are removed; different coding models
    at the same locus are retained as alternatives with provisional confidence.
    We do not collapse uncertain fragment containment or contradictory categories.
    """
    ordered=sorted(rows,key=lambda r:(decisions[r['candidate_id']]['category']=='complete',decisions[r['candidate_id']]['confidence']=='reviewed_complete',
                  decisions[r['candidate_id']]['model']['has_start'],
                  len(decisions[r['candidate_id']]['model']['full_ref_support']),
                  decisions[r['candidate_id']]['model']['or_bits']),reverse=True)
    clusters=[]
    for row in ordered:
        cid=row['candidate_id']; d=decisions[cid]; m=d['model']
        match=None
        for cluster in clusters:
            anchor=cluster[0]; ar=anchor['row']; ad=anchor['decision']; am=ad['model']
            identical_kind=ad['category']==d['category'] and same_locus(am,m)
            contained=(allow_containment and ad['category']=='complete' and d['category']=='partial'
                  and review.coding_intact(m) and m['canonical'] and m['protein'] in am['protein']
                  and review.overlap_fraction(m['exons'],am['exons'])>=.98 and shared_reading_frame(am,m))
            terminal=(allow_containment and ad['category']=='complete' and d['category']=='partial' and terminal_equivalent(am,m))
            fragment=(allow_containment and ad['category']=='partial' and d['category']=='partial' and fragment_alias(am,m))
            complete_fragment=(allow_containment and ad['category']=='complete' and d['category']=='partial' and fragment_alias(am,m))
            if (ar['sequence_id']==row['sequence_id'] and am['strand']==m['strand'] and (identical_kind or contained or terminal or fragment or complete_fragment)):
                rule='same_category_cds_overlap' if identical_kind else 'contained_intact_fragment' if contained else 'same_frame_terminal_variant' if terminal else 'same_frame_near_contained_complete_fragment' if complete_fragment else 'same_frame_near_contained_fragment'
                match=cluster; break
        member=dict(row=row,decision=d,alias_rule=rule if match is not None else 'representative')
        if match is None: clusters.append([member])
        else: match.append(member)
    final=outdir/'nonredundant'; final.mkdir(exist_ok=True)
    alias=[]; catalogue=[]; gff=['##gff-version 3']
    for cluster in clusters:
        anchor=cluster[0]; r=anchor['row']; d=anchor['decision']; m=d['model']; cid=r['candidate_id']
        different=len({x['decision']['model']['cds'] for x in cluster if x['decision']['category']==d['category']})>1
        confidence=d['confidence']
        if different and d['category']=='complete': confidence='provisional_complete'
        catalogue.append(dict(candidate_id=cid,aliases=','.join(x['row']['candidate_id'] for x in cluster[1:]),
                   sequence_id=r['sequence_id'],strand=m['strand'],category=d['category'],confidence=confidence,
                   length_aa=m['length'],distinct_alternative_cds_models=int(different),
                   locus_status='alternative_CDS_models' if different else 'representative',
                   start=min(s for s,e in m['exons']),end=max(e for s,e in m['exons']),
                   has_start=int(m['has_start']),has_stop=int(m['has_stop']),internal_stops=len(m['stops']),frameshifts=m['frameshifts']))
        for x in cluster:
            alias.append(dict(input_candidate_id=x['row']['candidate_id'],representative_id=cid,
                         retained=int(x is anchor),category=d['category'],input_model_category=x['decision']['category'],different_coding_models=int(different),alias_rule=x['alias_rule']))
        attrs=f"ID={cid};category={d['category']};confidence={confidence}"
        gff.append('\t'.join(map(str,[r['sequence_id'],'OR_model_review','mRNA',catalogue[-1]['start'],catalogue[-1]['end'],'.',m['strand'],'.',attrs])))
        cumulative=0
        for i,(s,e) in enumerate(m['exons'],1):
            gff.append('\t'.join(map(str,[r['sequence_id'],'OR_model_review','CDS',s,e,'.',m['strand'],(3-cumulative%3)%3,f'ID={cid}.CDS{i};Parent={cid}'])))
            cumulative+=e-s+1
    review.write_table(final/'OR_catalogue.tsv',catalogue,list(catalogue[0]))
    review.write_table(final/'OR_candidate_aliases.tsv',alias,list(alias[0]))
    residual=[]
    for index,left in enumerate(catalogue):
        lm=decisions[left['candidate_id']]['model']
        for right in catalogue[index+1:]:
            if left['sequence_id']!=right['sequence_id'] or left['strand']!=right['strand']: continue
            rm=decisions[right['candidate_id']]['model']
            a=review.overlap_fraction(lm['exons'],rm['exons']); b=review.overlap_fraction(rm['exons'],lm['exons'])
            if min(a,b)>.5:
                residual.append(dict(candidate_a=left['candidate_id'],candidate_b=right['candidate_id'],
                     category_a=left['category'],category_b=right['category'],cds_overlap_a=round(a,4),cds_overlap_b=round(b,4),
                     disposition='unresolved_overlap_not_merged'))
    review.write_table(final/'OR_residual_overlap.tsv',residual,list(residual[0]) if residual else ['candidate_a','candidate_b'])
    conflicted={r[k] for r in residual for k in ['candidate_a','candidate_b']}
    for r in catalogue:
        if r['candidate_id'] in conflicted:
            r['locus_status']='unresolved_overlap'
            if r['category']=='complete':r['confidence']='provisional_complete'
    bycid={r['candidate_id']:r for r in catalogue}
    for i,line in enumerate(gff):
        f=line.split('\t')
        if len(f)==9 and f[2]=='mRNA':
            cid=f[8].split(';')[0].removeprefix('ID=')
            f[8]=f"ID={cid};category={bycid[cid]['category']};confidence={bycid[cid]['confidence']};locus_status={bycid[cid]['locus_status']}"
            gff[i]='\t'.join(f)
    # Re-write after residual-conflict downgrades so TSV/GFF/FASTA agree.
    review.write_table(final/'OR_catalogue.tsv',catalogue,list(catalogue[0]))
    (final/'OR_models.gff3').write_text('\n'.join(gff)+'\n')
    for kind in ['complete','partial','pseudogene']:
        members=[r for r in catalogue if r['category']==kind]
        review.write_table(final/f'OR_{kind}.tsv',members,list(catalogue[0]))
        for field,extension in [('protein','proteins.faa'),('cds','cds.fna')]:
            write_fasta(final/f'OR_{kind}.{extension}',[FastaRecord(r['candidate_id'],f"category={kind} confidence={r['confidence']}",decisions[r['candidate_id']]['model'][field]) for r in members])
    for field,extension in [('protein','proteins.faa'),('cds','cds.fna')]:
        reviewed=[r for r in catalogue if r['confidence']=='reviewed_complete']
        write_fasta(final/f'OR_complete_reviewed.{extension}',[FastaRecord(r['candidate_id'],'',decisions[r['candidate_id']]['model'][field]) for r in reviewed])
    summary=dict(input_records=len(rows),representative_loci=len(clusters),representative_model_records=len(clusters),
          records_are_not_confirmed_gene_counts=True,aliased_records=len(rows)-len(clusters),
          category_counts=dict(Counter(r['category'] for r in catalogue)),confidence_counts=dict(Counter(r['confidence'] for r in catalogue)),
          alternative_model_loci=sum(r['distinct_alternative_cds_models'] for r in catalogue),
          residual_overlapping_pairs=len(residual),
          policy='Same contig/strand/category; reciprocal CDS overlap >=0.8 or >=0.65 with shared splice structure; quality-first non-transitive anchors; alternatives retained in final/; unresolved overlaps downgraded'+('; intact same-frame partial CDS >=98% contained AND peptide contained, OR intact OR same-frame/splice terminal variant with reciprocal CDS >=98% and endpoints within 12 nt; intact compatible partials may alias with smaller CDS >=90% contained, larger >=55% shared and <=30 nt outside; no damage/fusion/repeated-domain signal for new fragment/terminal rules' if allow_containment else ''),
          limitations='Other partial-overlap/contained fragment records remain unresolved; representative loci are not experimentally validated gene counts')
    (final/'OR_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


