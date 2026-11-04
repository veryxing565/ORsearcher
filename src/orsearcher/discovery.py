from __future__ import annotations
import concurrent.futures as futures,json
from collections import defaultdict
from pathlib import Path
from .fasta import read_fasta,FastaRecord,write_fasta
from .models import parse_miniprot_gff,attach_reference_lengths
from .support import review_or_catalogue as review
from .runtime import execute,save,LOG,digest

def cds_intervals(c):
    return [(f.start+c._offset,f.end+c._offset) for f in c.features if f.feature_type=='CDS']

def pool_alignments(raw):
    """Coding-overlap anchors, retaining members; never joins anchors transitively."""
    buckets=defaultdict(list);groups=[]
    for c in sorted(raw,key=lambda c:c.quality_key(),reverse=True):
        exons=cds_intervals(c)
        if not exons:continue
        bins=range(c.original_start//200000,c.original_end//200000+1);indices=set()
        for b in bins:indices.update(buckets[(c.original_id,c.strand,b)])
        match=None
        for i in sorted(indices):
            a=groups[i][0]
            if min(review.overlap_fraction(exons,cds_intervals(a)),review.overlap_fraction(cds_intervals(a),exons))>=.7:
                match=i;break
        if match is None:
            match=len(groups);groups.append([c])
            for b in bins:buckets[(c.original_id,c.strand,b)].append(match)
        else:groups[match].append(c)
    return groups

def align_pass(args,shards,references,folder,second=False):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    from .tuning import policy
    t=policy(args=args)
    def worker(shard):
        index=shard.with_suffix('.mpi')
        execute(['miniprot','-t',args.threads//args.jobs,'-d',index,shard],Path(str(index)+'.build.log'),inputs=[shard],resume=index.exists() and index.stat().st_size>0)
        if not index.is_file() or not index.stat().st_size:raise RuntimeError(f'No usable index {index}')
        gff=folder/(shard.stem+'.gff3')
        execute(['miniprot','-t',args.threads//args.jobs,'--gff','-G',args.max_intron,'--outc',t['second_output_coverage'] if second else t['first_output_coverage'],'--outs',t['alignment_output_score'],
                 '-N',args.max_hits,'--outn',args.max_hits,'-n',t['second_seed_hits'] if second else t['first_seed_hits'],'-p',t['second_alignment_ratio'] if second else t['first_alignment_ratio'],index,references],gff,inputs=[index,references])
        LOG.debug('%s finished %s',folder.name,shard.name);return gff
    with futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        jobs = {executor.submit(worker, shard): index for index, shard in enumerate(shards.shards)}
        results = [None] * len(jobs)
        for count, job in enumerate(futures.as_completed(jobs), 1):
            results[jobs[job]] = job.result()
            if count == 1 or count % 10 == 0 or count == len(jobs):
                LOG.info('Genome pass %d: %d/%d shards aligned', 2 if second else 1, count, len(jobs))
        return results

def read_pass(args,shards,gffs,refs,folder,label):
    checkpoint=Path(folder)/'models.json'
    if checkpoint.exists():return json.loads(checkpoint.read_text())
    mapping={s.virtual_id:s for s in shards.segments};raw=[]
    for gff in gffs:
        records=parse_miniprot_gff(gff,mapping);attach_reference_lengths(records,{k:len(v.sequence) for k,v in refs.items()})
        for c in records:c._offset=mapping[c.virtual_id].original_start
        raw.extend(records)
    LOG.info('%s parsing: raw=%d',label,len(raw));groups=pool_alignments(raw);byshard=defaultdict(list)
    for g in groups:byshard[g[0].shard].append(g)
    models=[]
    for shard in shards.shards:
        if shard.name not in byshard:continue
        sequences={r.id:r.sequence for r in read_fasta(shard)}
        for group in byshard[shard.name]:
            c=group[0];seq=sequences[c.virtual_id];features=sorted([f for f in c.features if f.feature_type=='CDS'],key=lambda f:f.start,reverse=c.strand=='-')
            cid='ORS2_'+args.sample+'_'+digest([c.original_id,c.strand,cds_intervals(c)])[:12]
            m=review.build_model(cid,'original',seq,c._offset,[(f.start,f.end) for f in features],c.strand,
                  int(features[0].phase) if features[0].phase.isdigit() else 0,c.frameshifts,c.reference_id,c.reference_start,c.reference_end,c.reference_length,c.identity)
            m.update(contig=c.original_id,protein_match=True,overlap=1.,search_passes=[label],
               search_references=sorted({x.reference_id for x in group}),alignment=dict(qstart=c.reference_start,qend=c.reference_end,qlen=c.reference_length),
               raw_alignment_count=len(group),annotation_gene_ids=[])
            variants=review.terminal_variants(m,seq,c._offset,refs)
            intact=[v for v in variants if review.coding_intact(v) and v['has_start'] and v['has_stop'] and v['canonical']]
            if intact:
                best=min(intact,key=lambda v:abs(v.get('terminal',{}).get('n_difference',999))+abs(v.get('terminal',{}).get('c_difference',999)))
                best.update({k:m[k] for k in ['contig','protein_match','search_passes','search_references','alignment','raw_alignment_count','annotation_gene_ids']})
                best['overlap']=1.;best['source']='original';m=best
            models.append(m)
        LOG.debug('%s extracted %s models=%d',label,shard.name,len(models))
    save(checkpoint,models);save(Path(folder)/'summary.json',dict(raw_alignments=len(raw),pooled_search_records=len(models),label=label))
    return models

def union_models(first,second):
    """Keep first-only models; collapse exact technical duplicates, not loci."""
    combined={}
    for m in [*first,*second]:
        key=(m['contig'],m['strand'],tuple(map(tuple,m['exons'])),m['cds'],m['frameshifts'])
        if key in combined:
            old=combined[key];old['search_passes']=sorted(set(old['search_passes']+m['search_passes']))
            old['search_references']=sorted(set(old['search_references']+m['search_references']))
            old.setdefault('discovery_alternatives',[]).extend(m.get('discovery_alternatives',[]))
            # Identical coordinates from different aligners/references remain
            # separate evidence, even when they share a technical representative.
            if (old.get('ref'),old.get('alignment_tool','miniprot'))!=(m.get('ref'),m.get('alignment_tool','miniprot')):
                old['discovery_alternatives'].append({k:v for k,v in m.items() if k!='discovery_alternatives'})
        else:combined[key]=m
    return list(combined.values())

def adaptive_seeds(models,args=None):
    from .tuning import policy
    t=policy(models[0] if models else None,args=args)
    return [m for m in models if review.coding_intact(m) and m['has_start'] and m['has_stop'] and m['canonical']
            and m.get('family')=='OR' and (m['or_ga'] or m['or_evalue']<=1e-3 and m['or_bits']>=100)
            and (bool(m.get('full_ref_support')) or not t['adaptive_require_full_reference']) and not m.get('repeated_or_domains',False)
            and m['ref_cov']>=t['adaptive_reference_coverage'] and m['identity']>=t['adaptive_min_identity']/100 and t['adaptive_min_aa']<=m['length']<=t['adaptive_max_aa']]
