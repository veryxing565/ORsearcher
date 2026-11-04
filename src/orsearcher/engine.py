from __future__ import annotations
import concurrent.futures as futures,csv,json,logging,os,re,shutil,subprocess,time
from collections import Counter,defaultdict
from pathlib import Path
from types import SimpleNamespace
from . import __version__
from .fasta import FastaRecord,read_fasta,write_fasta,translate,reverse_complement
from .sharding import create_shards
from .models import parse_attributes
from .discovery import align_pass,read_pass,union_models,adaptive_seeds
from .translated import discover as translated_discover
from .recall import discover as recall_discover,export_candidates
from .evidence import Evidence,route
from .references import audit_references
from .model_store import save_groups
from .runtime import LOG,digest,metadata,save,execute,fasta_if_any,file_digest
from .support import review_or_catalogue as review,audit_or_loci as audit,crosscheck_or_catalogue as cross
from .support.finalize_catalogue_review import export_unique

def table(path):
    with Path(path).open(encoding='utf-8') as handle:return list(csv.DictReader(handle,delimiter='\t'))

def load_samples(args):
    """Resolve manifest-relative inputs against the manifest, not the CWD."""
    if not args.manifest:
        return [dict(sample=args.sample,genome=str(args.genome.resolve()),
                     annotation=str(args.annotation.resolve()) if args.annotation else '')]
    rows=table(args.manifest)
    root=args.manifest.resolve().parent
    for row in rows:
        for name in ['genome','annotation']:
            if row.get(name):
                path=Path(row[name]).expanduser()
                row[name]=str((root/path if not path.is_absolute() else path).resolve())
    return rows

def discovery_reports(out):
    reports={}
    warnings=[]
    for stage in ['translated','recall']:
        path=Path(out)/'discovery'/stage/'summary.json'
        if path.is_file():
            reports[stage]=json.loads(path.read_text(encoding='utf-8'))
            warnings.extend(stage+':'+item for item in reports[stage].get('warnings',[]))
            if reports[stage].get('status') in ['error','failed'] and not reports[stage].get('warnings'):
                warnings.append(stage+':stage_failed')
    return reports,warnings

def row_for(m):
    return dict(candidate_id=m['cid'],sequence_id=m['contig'],strand=m['strand'],start=min(s for s,e in m['exons']),end=max(e for s,e in m['exons']),
       protein_length=m['length'],frameshifts=m['frameshifts'],original_category='pseudogene' if m['frameshifts'] or m['stops'] else 'complete' if m['has_start'] and m['has_stop'] else 'partial')

def dependency_check(args):
    polishing=getattr(args,'command',None)=='polish'
    required=['miniprot','diamond','hmmsearch','samtools','exonerate']
    if not polishing:
        required.extend(['gffread','mafft','trimal','FastTree','Rscript'])
        if getattr(args,'translated_search','auto')!='off':required.extend(['tblastn','makeblastdb'])
        if getattr(args,'recall_search','auto')!='off':required.extend(['metaeuk','genewise'])
    missing=[name for name in required if not shutil.which(name)]
    if missing:raise RuntimeError('Missing dependencies: '+', '.join(missing)+'; install environment.yml with mamba')
    versions={}
    for name in required:
        binary=shutil.which(name)
        if not binary:raise RuntimeError('Missing dependency: '+name)
        cmd=[binary,'version'] if name=='metaeuk' else [binary,'-version'] if name in ['tblastn','makeblastdb','genewise'] else [binary,'--version'] if name not in ['diamond','hmmsearch','samtools','FastTree','trimal'] else [binary,'version'] if name=='diamond' else [binary,'-h'] if name=='hmmsearch' else [binary,'--version'] if name=='samtools' else [binary,'-help'] if name=='FastTree' else [binary,'--version']
        try:
            p=subprocess.run(cmd,capture_output=True,text=True,errors='replace',timeout=10)
            if p.returncode<0:raise RuntimeError('Dependency '+name+' crashed during version check (signal '+str(-p.returncode)+'); check CPU/platform compatibility')
            versions[name]=dict(binary=binary,returncode=p.returncode,version=(p.stdout+p.stderr).splitlines()[:4])
        except subprocess.TimeoutExpired:versions[name]=dict(binary=binary,version_check='timeout')
    if not polishing:
        p=subprocess.run(['Rscript','-e','stopifnot(requireNamespace("ape",quietly=TRUE))'],capture_output=True,text=True,timeout=30)
        if p.returncode:raise RuntimeError('R package ape is required: '+p.stderr)
    if not polishing and getattr(args,'recall_search','auto')!='off':
        cfg=getattr(args,'wise_config',None) or Path(versions['genewise']['binary']).resolve().parents[1]/'share/wise2/wisecfg'
        if not (cfg/'codon.table').is_file() or not (cfg/'BLOSUM45.bla').is_file():raise RuntimeError('Missing GeneWise config files: '+str(cfg))
        versions['genewise']['config_directory']=str(cfg)
    return versions

def prepare_references(args,root):
    folder=root/'shared/reference';folder.mkdir(parents=True,exist_ok=True)
    refs={};excluded=[]
    complete_path=getattr(args,'complete_reference',None)
    complete={r.id:r for r in read_fasta(complete_path)} if complete_path else {}
    exclude_sequences={r.sequence.rstrip('*') for r in read_fasta(args.exclude_reference_proteins)} if args.exclude_reference_proteins else set()
    for path,family in [(args.or_reference,'or'),(args.gr_reference,'gr')]:
        records=list(read_fasta(path));seen=set();selected=[]
        if family=='or':
            originals={r.id:r.sequence.rstrip('*') for r in records}
            for rid,record in complete.items():
                if originals.get(rid)!=record.sequence.rstrip('*'):
                    raise ValueError('Complete reference must have matching ID AND sequence in OR reference: '+rid)
        for r in records:
            if r.id in seen:raise ValueError('Duplicate reference ID '+r.id)
            seen.add(r.id)
            if any(re.search(pattern,r.description) for pattern in args.exclude_reference_pattern) or r.sequence.rstrip('*') in exclude_sequences:
                excluded.append(dict(reference_id=r.id,family=family,reason='held_out_header_or_exact_sequence'));continue
            if len(r.sequence)<50 or '*' in r.sequence.rstrip('*'):raise ValueError('Invalid reference protein '+r.id)
            selected.append(FastaRecord(r.id,r.description,r.sequence.rstrip('*')))
        if not selected:raise ValueError('Reference exclusions left no '+family+' proteins')
        fasta=folder/(family+'.faa')
        if not fasta.exists():write_fasta(fasta,selected)
        db=folder/family
        execute(['diamond','makedb','--in',fasta,'--db',db,'--threads',args.threads],folder/(family+'.makedb.log'),inputs=[fasta],resume=Path(str(db)+'.dmnd').exists())
        setattr(args,family+'_db',db)
        if family=='or':
            refs={r.id:r for r in selected}
            audit_references(args,folder,selected)
            audited = {rid for rid, meta in args.reference_metadata.items() if meta.get('reconstruction_allowed')}
            args.complete_reference_ids = set(complete) & set(refs) if complete_path else audited
            if complete_path and getattr(args, 'use_added_templates', False):
                args.complete_reference_ids |= audited & getattr(args, 'user_or_ids', set())
            LOG.info('OR references: %d protein(s), %d eligible completeness template(s)',
                     len(refs), len(args.complete_reference_ids))
        else:
            LOG.info('GR references: %d protein(s)', len(selected))
    save(folder/'excluded_references.json',excluded)
    return refs

def annotation_references(args,out,evidence):
    if args.annotation is None:return []
    folder=out/'annotation';folder.mkdir(exist_ok=True)
    execute(['gffread',args.annotation,'-g',args.genome,'-y',folder/'proteins.faa','-x',folder/'cds.fna'],folder/'extract.log',inputs=[args.annotation,args.genome])
    models=[]
    for r in fasta_if_any(folder/'proteins.faa',read_fasta):
        if len(r.sequence)<100:continue
        models.append(dict(id=r.id,protein=r.sequence.rstrip('*'),length=len(r.sequence.rstrip('*'))))
    evidence.apply(models,'annotation')
    return [FastaRecord('annotation_'+str(i),'source='+m['id'],m['protein']) for i,m in enumerate(models) if route(m)=='OR_review' and m['family']=='OR']

def choose_refs(model,refs,limit=10,metadata=None,policy='audited'):
    names=list(model.get('full_ref_support',[]))
    hits=sorted(model.get('homology_hits',[]),key=lambda h:((h['se']-h['ss']+1)/h['slen']>=.8,h['bits']),reverse=True)
    names.extend(h['ref'] for h in hits if h['ref'] in refs and 100<=len(refs[h['ref']].sequence)<=1500)
    names.extend(r for r in model.get('search_references',[]) if r in refs)
    names.extend(a.get('ref','') for a in model.get('discovery_alternatives',[]))
    if model.get('ref') in refs:names.append(model['ref'])
    names=[r for r in dict.fromkeys(names) if r in refs and 100<=len(refs[r].sequence)<=1500
           and (metadata is None or metadata.get(r,{}).get('reconstruction_allowed',False) or policy=='active-external')]
    selected={};seen=set()
    for r in names:
        if refs[r].sequence in seen:continue
        selected[r]=refs[r];seen.add(refs[r].sequence)
        if len(selected)>=limit:break
    return selected

def initial_case(row,m,refs,run,args,lengths,rebuild):
    case=run/'cases'/row['candidate_id'];case.mkdir(parents=True,exist_ok=True)
    checkpoint=case/'models.json'
    status_path=case/'local_status.json'
    if checkpoint.exists() and status_path.exists() and json.loads(status_path.read_text()).get('status')=='completed':return json.loads(checkpoint.read_text())
    start=max(1,int(row['start'])-args.flank_bp);end=min(lengths[row['sequence_id']],int(row['end'])+args.flank_bp)
    if rebuild and refs:
        span=int(row['end'])-int(row['start'])+1
        if span>args.max_window_bp:raise ValueError('Structural model span exceeds extraction limit; candidate retained unresolved')
        bounded_flank=min(args.flank_bp,max(0,(args.max_window_bp-span)//2))
        settings=dict(vars(args));settings.update(flank_bp=bounded_flank,threads_per_job=1,samtools='samtools',miniprot='miniprot',outdir=run,miniprot_timeout=args.local_timeout)
        params=SimpleNamespace(**settings)
        values=review.model_candidates(row,case,[(s,e,0) for s,e in m['exons']],m['protein'],refs,lengths,params,original_model=m)
    else:
        result=subprocess.run(['samtools','faidx',str(args.genome),f"{row['sequence_id']}:{start}-{end}"],capture_output=True,text=True,check=True)
        seq=''.join(result.stdout.splitlines()[1:]).upper();offset=start-1
        write_fasta(case/'region.fna',[FastaRecord('region','',seq)]);write_fasta(case/'references.faa',refs.values())
        save(case/'region.json',dict(offset=offset,end=end,contig=row['sequence_id'],strand=row['strand']))
        (case/'miniprot.gff3').write_text('##gff-version 3\n')
        original=dict(m,source='original',protein_match=True,overlap=1.);values=[original,*review.terminal_variants(original,seq,offset,refs)]
    # Discovery pooling is not a structural verdict. Carry distinct reference,
    # tool and splice alternatives forward instead of silently discarding them.
    original=next(v for v in values if v['source']=='original')
    seen={(tuple(map(tuple,v['exons'])),v.get('ref',''),audit.tool(v)) for v in values}
    for alt in m.get('discovery_alternatives',[]):
        signature=(tuple(map(tuple,alt['exons'])),alt.get('ref',''),audit.tool(alt))
        if signature in seen:continue
        v={k:value for k,value in alt.items() if k!='discovery_alternatives'}
        v.update(cid=m['cid'],id=m['cid']+'.discovery_'+digest([signature,alt['cds']])[:12],source='discovery_'+alt['source'],
                 overlap=review.overlap_fraction(original['exons'],alt['exons']))
        values.append(v);seen.add(signature)
    for v in values:v.update(contig=row['sequence_id'],annotation_gene_ids=[],fusion_pairs=[])
    save(checkpoint,values);save(case/'local_status.json',dict(status='completed',local_rebuilt=bool(rebuild and refs),reference_count=len(refs)))
    return values


def initialize_fusion_evidence(models):
    """New aligner models need spatial flags BEFORE early rescue decisions."""
    anchors=audit.anchor_groups(models)
    for m in models:m['fusion_pairs']=audit.fusion_links(m,anchors)
    return anchors

def structure(args,models,out,evidence,refs):
    folder=out/'structure';folder.mkdir(exist_ok=True);rows=[row_for(m) for m in models]
    if len({r['candidate_id'] for r in rows})!=len(rows):raise ValueError('Candidate ID collision')
    if not Path(str(args.genome)+'.fai').exists():
        execute(['samtools','faidx',args.genome],out/'genome_index.log',inputs=[args.genome],resume=False)
    lengths=review.read_fai(Path(str(args.genome)+'.fai'))
    global_anchors=audit.anchor_groups(models)
    planned={}
    for row,m in zip(rows,models):
        fusion=audit.fusion_links(m,global_anchors)
        robust=(review.coding_intact(m) and m['has_start'] and m['has_stop'] and m['canonical'] and len(m['full_ref_support'])>=2
               and m['or_ga'] and 250<=m['length']<=650 and not fusion and not m['repeated_or_domains'])
        planned[m['cid']]=dict(model=m,references=choose_refs(m,refs,limit=getattr(args,'structure_references',10),metadata=args.reference_metadata,policy=getattr(args,'reconstruction_reference_policy','audited')),rebuild=not robust or getattr(args,'force_structural_review',False))
    byid={};failures=[]
    with futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(initial_case,r,planned[r['candidate_id']]['model'],planned[r['candidate_id']]['references'],folder,args,lengths,planned[r['candidate_id']]['rebuild']):r for r in rows}
        for i,job in enumerate(futures.as_completed(jobs),1):
            r=jobs[job];cid=r['candidate_id']
            try:byid[cid]=job.result()
            except Exception as error:
                LOG.exception('Local case failed %s',cid);failures.append(dict(candidate_id=cid,status='error',error=str(error)))
                m=dict({k:v for k,v in planned[cid]['model'].items() if k!='discovery_alternatives'},source='original',protein_match=True,fusion_pairs=[]);byid[cid]=[m]
            if i%20==0 or i==len(rows):LOG.info('Local structures %d/%d',i,len(rows))
    initial=[m for values in byid.values() for m in values]
    neighbors=audit.collect_neighbors([r for r in rows if (folder/'cases'/r['candidate_id']/'region.fna').exists()],folder,refs)
    evidence.apply(initial+neighbors,'local_models')
    second=folder/'crosscheck';second.mkdir(exist_ok=True);statuses=[];extras=[]
    parameters=SimpleNamespace(outdir=second,timeout=args.exonerate_timeout,exonerate='exonerate',max_intron=args.max_intron,max_window_bp=args.max_window_bp,references_per_candidate=3,pseudogene_references=4,exonerate_word_length=args.exonerate_word_length,legacy_rescue=args.exonerate_legacy_rescue)
    from .structure_rescue import legacy_priority
    legacy_allowed={r['candidate_id'] for r in sorted(rows,
        key=lambda r:(legacy_priority(byid[r['candidate_id']]),r['candidate_id']),reverse=True)[:args.exonerate_legacy_budget]}
    def second_job(row):
        cid=row['candidate_id'];case=second/'cases'/cid
        if (case/'status.json').exists() and json.loads((case/'status.json').read_text()).get('status')=='completed' and all((case/p).exists() for p in ['models.json','miniprot_frame_events.json','exonerate_frame_events.json']):
            return json.loads((case/'models.json').read_text()),json.loads((case/'status.json').read_text())
        if not (folder/'cases'/cid/'region.fna').exists():return [],dict(candidate_id=cid,status='local_failed')
        chosen=max(byid[cid],key=review.model_rank)
        local_params=SimpleNamespace(**vars(parameters),allow_legacy=cid in legacy_allowed)
        return cross.process(row,byid[cid],folder,second,dict(selected_model_id=chosen['id'],model_category=row['original_category']),local_params)
    with futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(second_job,r):r for r in rows}
        for i,job in enumerate(futures.as_completed(jobs),1):
            r=jobs[job]
            try:extra,status=job.result()
            except Exception as error:LOG.exception('Second aligner failed');extra=[];status=dict(candidate_id=r['candidate_id'],status='error',error=str(error))
            extras.extend(extra);statuses.append(status)
            if status['status']!='completed':LOG.debug('Second aligner case: %s',status)
            if i%100==0 or i==len(rows):LOG.info('Second aligner %d/%d',i,len(rows))
    evidence.apply(extras,'exonerate_models')
    for m in extras:m.update(contig=next(r['sequence_id'] for r in rows if r['candidate_id']==m['cid']),annotation_gene_ids=[])
    for m in extras:byid[m['cid']].append(m)
    if getattr(args,'identity_policy','legacy')=='multichannel':
        from .structure_rescue import reconstruct
        initialize_fusion_evidence([m for values in byid.values() for m in values]+neighbors)
        reconstruct(args,rows,byid,folder,evidence,refs)
    allmodels=[m for v in byid.values() for m in v]+neighbors
    for m in allmodels:m.setdefault('annotation_gene_ids',[])
    anchors=audit.anchor_groups(allmodels)
    for m in allmodels:m['fusion_pairs']=audit.fusion_links(m,anchors)
    # Multi-HSP repeat analysis is bounded to flagged models, not all proteins.
    flagged=[m for m in allmodels if m['fusion_pairs'] or m['repeated_or_domains'] or m['length']>650]
    if flagged:
        repeat_args=SimpleNamespace(outdir=folder,workers=args.workers,or_hmm=args.or_hmm,or_db=args.or_db,diamond='diamond',hmmsearch='hmmsearch')
        cache=folder/'protein_structure/flags.json'
        if cache.exists():
            repeat_flags=json.loads(cache.read_text())
            for m in flagged:
                if m['id'] in repeat_flags:m.update(repeat_flags[m['id']])
        else:
            audit.protein_repeats(flagged,repeat_args)
            save(cache,{m['id']:{k:m[k] for k in ['repeated_or_domains','repeated_reference_regions']} for m in flagged})
    if args.annotation:audit.annotation_context(args.annotation,allmodels)
    # Share alternatives only when substantial coding/genomic overlap exists.
    index=defaultdict(list)
    for m in allmodels:
        if m['source']!='original':index[(m['contig'],m['strand'])].append(m)
    for row in rows:
        cid=row['candidate_id'];original=next(m for m in byid[cid] if m['source']=='original')
        seen={(tuple(map(tuple,m['exons'])),m['ref'],audit.tool(m)) for m in byid[cid]}
        for m in index[(row['sequence_id'],row['strand'])]:
            signature=(tuple(map(tuple,m['exons'])),m['ref'],audit.tool(m))
            if signature in seen:continue
            n=dict(m,overlap=review.overlap_fraction(original['exons'],m['exons']))
            if audit.near_target(n,original):
                n.update(cid=cid,id=cid+'.shared_'+digest(m['id'])[:12],source='transferred_'+m['source']);byid[cid].append(n);seen.add(signature)
    frames=defaultdict(list)
    for row in rows:
        case=second/'cases'/row['candidate_id']
        for name,aligner in [('miniprot_frame_events.json','miniprot'),('exonerate_frame_events.json','exonerate')]:
            if (case/name).exists():frames[row['sequence_id']].append((aligner,json.loads((case/name).read_text())))
        rescue_frames=folder/'rescue/cases'/row['candidate_id']/'exonerate_frame_events.json'
        if rescue_frames.exists():frames[row['sequence_id']].append(('exonerate',json.loads(rescue_frames.read_text())))
    decisions={};damage_rows=[];qc=[];family_unresolved=[];accepted=[];selection_rows=[];partial_rows=[]
    for row in rows:
        cid=row['candidate_id'];original=next(m for m in byid[cid] if m['source']=='original')
        damage=audit.damage_summary(cid,byid[cid],original,second,frames[row['sequence_id']])
        d=audit.decide_locus(row,byid[cid],damage);m=d['model']
        if d['category']=='partial':
            viable=[x for x in byid[cid] if audit.near_target(x,original) and audit.structure_candidate(x)
                    and x['canonical'] and not x['fusion_pairs']
                    and not(x.get('repeated_or_domains') and x.get('repeated_reference_regions'))]
            for option in viable or [m]:
                partial_rows.append(dict(candidate_id=cid,model_id=option['id'],selected=int(option['id']==m['id']),
                    source=option['source'],length_aa=option['length'],canonical=option['canonical'],
                    coding_intact=review.coding_intact(option),rank=str(audit.partial_selection_rank(option,viable)),
                    **audit.partial_selection_evidence(option,viable)))
        if m['family']!='OR':
            family_unresolved.append(dict(candidate_id=cid,family=m['family'],or_bits=m['or_bits'],gr_bits=m['gr_bits'],reason=d['reason']))
            if getattr(args,'identity_policy','legacy')!='multichannel':continue
        accepted.append(row);decisions[cid]=d
        good=[x for x in byid[cid] if audit.same_region(x,original) and not x['fusion_pairs']
              and not(x['repeated_or_domains'] and x['repeated_reference_regions']) and review.complete_model(dict(x,overlap=1.,bridge=False))]
        support=audit.compatible_support(m,byid[cid],independent=True)
        selection_rows.append(dict(candidate_id=cid,selected_model_id=m['id'],category=d['category'],reference=m.get('ref',''),
             reference_origin=m.get('alignment_reference_origin','unknown'),reference_qualified=int(audit.external_alignment(m)),
             true_alignment_tool=audit.tool(m),external_alignment_identity=audit.external_identity(m),reference_coverage=m.get('ref_cov',0),
             qualified_exact_protein_match=int(audit.external_exact(m)),independent_tool_reference_sequence_pairs=len(support),
             complete_selection_rank=str(audit.selection_rank(m,good)) if good else 'partial_fallback',reason=d['reason']))
        selection_rows[-1]['partial_selection_rank']=str(audit.partial_selection_rank(m,
            [x for x in byid[cid] if audit.near_target(x,original) and audit.structure_candidate(x) and x['canonical']
             and not x['fusion_pairs']])) if d['category']=='partial' else ''
        if damage:damage_rows.append(dict(candidate_id=cid,kind=damage['kind'],coordinates=str(damage['coordinates']),references=','.join(damage['references']),tools=','.join(damage['tools']),category=d['category']))
        qc.append(dict(candidate_id=cid,category=d['category'],confidence=d['confidence'],has_start=m['has_start'],has_stop=m['has_stop'],internal_stops=len(m['stops']),frameshifts=m['frameshifts'],canonical=m['canonical'],
                 fusion_pairs=str(m['fusion_pairs']),local_rebuilt=planned[cid]['rebuild'],original_search_passes=','.join(planned[cid]['model']['search_passes']),reason=d['reason']))
    for name,values,header in [('integrity_qc.tsv',qc,['candidate_id']),('model_selection.tsv',selection_rows,['candidate_id']),
                            ('partial_selection_audit.tsv',partial_rows,['candidate_id','model_id','selected']),
                            ('disruption_evidence.tsv',damage_rows,['candidate_id']),('family_unresolved.tsv',family_unresolved,['candidate_id']),('local_failures.tsv',failures,['candidate_id'])]:
        review.write_table(folder/name,values,list(values[0]) if values else header)
    save_groups(folder/'all_models.json',byid);save(folder/'crosscheck_summary.json',dict(statuses=statuses,counts=dict(Counter(s['status'] for s in statuses)),legacy_status_counts=dict(Counter(s.get('legacy_status','not_recorded') for s in statuses)),
             candidate_count=len(rows),local_rebuilt=sum(x['rebuild'] for x in planned.values()),neighbor_models=len(neighbors),spatial_anchors=len(anchors)))
    return accepted,byid,decisions

def phylogeny(args,rows,decisions,refs,out):
    folder=out/'phylogeny';folder.mkdir(exist_ok=True);targets={r['candidate_id']:decisions[r['candidate_id']]['model'] for r in rows}
    target_sequences={m['protein'] for m in targets.values()};records=[];labels=[];excluded=[]
    # Honor held-out exclusions in diagnostic references as well as discovery.
    for family,path in [('OR',Path(str(args.or_db)+'.faa')),('GR',Path(str(args.gr_db)+'.faa'))]:
        for i,r in enumerate(read_fasta(path)):
            tid='ref'+family+'-'+''.join(c if c.isalnum() or c in '_.|-' else '_' for c in r.id)
            if r.sequence in target_sequences:excluded.append(dict(reference_id=r.id,family=family,reason='identical_to_target'));continue
            records.append(FastaRecord(tid,'',r.sequence));labels.append(dict(tree_id=tid,original_id=r.id,family=family,state='reference'))
    omitted=[]
    for cid,m in targets.items():
        if m['length']<150:omitted.append(dict(candidate_id=cid,reason='shorter_than_150aa'));continue
        if m['length']>getattr(args,'phylogeny_max_aa',1500):omitted.append(dict(candidate_id=cid,reason='long_model_alignment_deferred'));continue
        records.append(FastaRecord(cid,'',m['protein'].replace('*','X')));labels.append(dict(tree_id=cid,original_id=cid,family='OR',state='candidate'))
    review.write_table(folder/'id_map.tsv',labels,list(labels[0]));review.write_table(folder/'excluded_identical_references.tsv',excluded,list(excluded[0]) if excluded else ['reference_id'])
    review.write_table(folder/'omitted_short.tsv',omitted,list(omitted[0]) if omitted else ['candidate_id'])
    write_fasta(folder/'input.faa',records)
    if not targets:return {}
    timeout=getattr(args,'phylogeny_timeout',900)
    execute(['mafft','--auto','--thread',args.threads,folder/'input.faa'],folder/'aligned.faa',inputs=[folder/'input.faa'],timeout=timeout,kill_process_group=True)
    execute(['trimal','-in',folder/'aligned.faa','-out',folder/'trimmed.faa','-gt','.2'],folder/'trim.log',inputs=[folder/'aligned.faa'],timeout=timeout)
    execute(['FastTree','-lg','-gamma',folder/'trimmed.faa'],folder/'tree.nwk',inputs=[folder/'trimmed.faa'],timeout=timeout)
    script=Path(__file__).parent/'tree_qc.R'
    execute(['Rscript',script,folder/'tree.nwk',folder/'id_map.tsv',folder/'trimmed.faa',folder],folder/'tree_qc.log',inputs=[script,folder/'tree.nwk',folder/'id_map.tsv',folder/'trimmed.faa'],timeout=timeout)
    flags={r['candidate_id']:r for r in table(folder/'candidate_qc.tsv')}
    for cid,r in flags.items():
        # Tree side/gaps alone do not delete or relabel an OR gene.
        if r['tree_side']=='GR_side' or float(r['gap_fraction'])>.5:
            if decisions[cid]['category']=='complete':decisions[cid]['confidence']='provisional_complete'
            decisions[cid]['reason']+=';phylogeny_review_flag'
    return flags

def validate(genome,folder):
    candidate=(folder/'candidate_catalogue.tsv').exists()
    catalogue=table(folder/('candidate_catalogue.tsv' if candidate else 'OR_catalogue.tsv'));proteins={};cdss={};features=defaultdict(list)
    if candidate:
        proteins={r.id:r.sequence for r in read_fasta(folder/'candidate_proteins.faa')};cdss={r.id:r.sequence for r in read_fasta(folder/'candidate_cds.fna')}
    for kind in ['complete','partial','pseudogene']:
        proteins.update({r.id:r.sequence for r in fasta_if_any(folder/f'OR_{kind}.proteins.faa',read_fasta)})
        cdss.update({r.id:r.sequence for r in fasta_if_any(folder/f'OR_{kind}.cds.fna',read_fasta)})
    for line in (folder/('candidate_models.gff3' if candidate else 'OR_models.gff3')).read_text().splitlines():
        if line.startswith('#'):continue
        f=line.split('\t')
        if f[2]=='CDS':features[parse_attributes(f[8])['Parent']].append((int(f[3]),int(f[4]),int(f[7])))
    from .indexed_genome import IndexedGenome
    with IndexedGenome(genome) as reader:
        for index,r in enumerate(catalogue,1):
            cid=r['candidate_id'];strand=r['strand'];ordered=sorted(features[cid],reverse=strand=='-');left=min(s for s,e,p in ordered);right=max(e for s,e,p in ordered)
            seq=reader.fetch(r['sequence_id'],left,right)
            parts=[];n=0
            for s,e,phase in ordered:
                if phase!=(3-n%3)%3:raise ValueError('Export phase mismatch '+cid)
                piece=seq[s-left:e-left+1];parts.append(reverse_complement(piece) if strand=='-' else piece);n+=e-s+1
            cds=''.join(parts);protein=translate(cds).removesuffix('*')
            if cds!=cdss[cid] or protein!=proteins[cid]:raise ValueError('Genome/GFF/CDS/protein mismatch '+cid)
            if r.get('category',r.get('model_category'))=='complete' and (not cds.startswith('ATG') or not translate(cds).endswith('*') or '*' in protein or 'X' in protein or len(cds)%3):raise ValueError('Complete-model failure '+cid)
            if index%10000==0 or index==len(catalogue):LOG.info('Export validation %d/%d (%s)',index,len(catalogue),folder.name)
    result=dict(checked_models=len(catalogue),coordinate_sequence_mismatches=0,complete_coding_failures=0)
    save(folder/'export_validation.json',result);return result

def run_sample(args,root,refs,evidence):
    start=time.monotonic();out=root/'samples'/args.sample;out.mkdir(parents=True,exist_ok=True)
    from .logging_utils import ConciseFormatter
    handler=logging.FileHandler(out/'run.log',encoding='utf-8');handler.setLevel(logging.INFO);handler.setFormatter(ConciseFormatter('%(asctime)s %(levelname)s %(message)s'));LOG.addHandler(handler)
    try:
        LOG.info('Sample %s: genome=%s annotation=%s',args.sample,args.genome,args.annotation)
        if not Path(str(args.genome)+'.fai').exists():execute(['samtools','faidx',args.genome],out/'genome_index.log',inputs=[args.genome],resume=False)
        if getattr(args,'command','run')=='polish':
            from .polish import polish
            result=polish(args,out,refs,evidence)
            summary=dict(sample=args.sample,version=__version__,status=result['status'],
                         seconds=round(time.monotonic()-start,2),analysis_scope='reviewed_models_consolidation_and_bounded_repair_not_fresh_genome_search',
                         polishing=result,output=str(out/'results'))
            save(out/'sample_summary.json',summary)
            return summary
        if getattr(args,'command','run')=='refine':
            from .refine import annotate
            expanded=[]
            for path in args.candidate_models:
                data=json.loads(path.read_text())
                if not isinstance(data,list):raise ValueError('Candidate JSON must contain a model list, not a GFF or decision dictionary')
                for m in data:expanded.extend([m,*m.get('discovery_alternatives',[])])
            evidence.apply(expanded,'cached_discovery_'+args.sample)
            pool=export_candidates(expanded,out/'recall_candidates');pool_validation=validate(args.genome,out/'recall_candidates')
            result=annotate(args,expanded,refs,evidence,out)
            summary=dict(sample=args.sample,version=__version__,status=result['status'],seconds=round(time.monotonic()-start,2),
                analysis_scope='cached_discovery_new_structure_identity_annotation_not_fresh_genome_search',candidate_pool=pool,candidate_validation=pool_validation,annotation=result,output=str(out))
            save(out/'sample_summary.json',summary);return summary
        LOG.info('[%s] Discovery: preparing overlap windows and packed shards',args.sample)
        shards=create_shards(args.genome,out/'shared_genome',args.chunk_size,args.overlap,args.shard_size,True)
        save(out/'genome_summary.json',dict(total_bases=shards.total_bases,input_sequences=shards.sequence_count,shards=len(shards.shards)))
        LOG.info('[%s] Genome: %s bases, %d sequence(s), %d shard(s)',
                 args.sample, f'{shards.total_bases:,}', shards.sequence_count, len(shards.shards))
        reference1=root/'shared/reference/or.discovery.faa'
        search_refs=dict(refs)
        extra_seeds={}
        if getattr(args,'search_seeds',None):
            for r in read_fasta(args.search_seeds):
                if len(r.sequence)<30 or '*' in r.sequence or not set(r.sequence)<=set('ACDEFGHIKLMNPQRSTVWYXBZ'):
                    raise ValueError('Invalid search-only seed '+r.id)
                sid='adaptive_external_'+digest(r.sequence)[:20]
                extra_seeds[sid]=FastaRecord(sid,'self_generated search_only=true',r.sequence)
            search_refs.update(extra_seeds)
            reference1=out/'discovery/first_search_reference.faa'
            write_fasta(reference1,[*args.discovery_references.values(),*extra_seeds.values()])
            save(out/'search_seed_audit.json',dict(input=str(args.search_seeds),seed_ids=sorted(extra_seeds),external_identity_reference_count=len(refs),seed_ids_in_identity_db=0,seed_ids_in_complete_templates=0))
        LOG.info('[%s] Discovery: first genome pass (%d shards)',args.sample,len(shards.shards))
        first_paths=align_pass(args,shards,reference1,out/'discovery/pass1')
        first=read_pass(args,shards,first_paths,search_refs,out/'discovery/pass1','first')
        evidence.apply(first,'first_'+args.sample)
        seeds=adaptive_seeds(first,args);annotation=annotation_references(args,out,evidence)
        ref2=dict(search_refs)
        for i,m in enumerate(seeds):ref2['adaptive_'+str(i)]=FastaRecord('adaptive_'+str(i),'self_generated',m['protein'])
        ref2.update({r.id:r for r in annotation})
        query_records=dict(args.discovery_references)
        query_records.update({k:v for k,v in ref2.items() if k not in refs})
        query2=out/'discovery/second_reference.faa';write_fasta(query2,query_records.values())
        LOG.info('Second pass: self_seeds=%d annotation=%d; zero seeds continues with original references',len(seeds),len(annotation))
        LOG.info('[%s] Discovery: relaxed second genome pass',args.sample)
        second_paths=align_pass(args,shards,query2,out/'discovery/pass2',True)
        second=read_pass(args,shards,second_paths,ref2,out/'discovery/pass2','second')
        union=union_models(first,second)
        evidence.apply(union,'pre_translated_'+args.sample)
        covered=[m for m in union if route(m)=='OR_review' and m['family']=='OR']
        LOG.info('[%s] Discovery: translated-search rescue',args.sample)
        translated=translated_discover(args,shards,refs,covered,out/'discovery/translated',evidence,reference1)
        union=union_models(union,translated)
        LOG.info('[%s] Discovery: MetaEuk/GeneWise rescue',args.sample)
        recall=recall_discover(args,shards,refs,union,out/'discovery/recall',evidence,out/'discovery/translated')
        union=union_models(union,recall)
        expanded=[x for m in union for x in [m,*m.get('discovery_alternatives',[])]]
        # Broad candidate export precedes family routing; no OR/GR deletion.
        pool=export_candidates(expanded,out/'recall_candidates')
        pool_validation=validate(args.genome,out/'recall_candidates')
        if args.recall_only:
            rs=json.loads((out/'discovery/recall/summary.json').read_text());ts=json.loads((out/'discovery/translated/summary.json').read_text())
            summary=dict(sample=args.sample,version=__version__,status='completed_with_warnings' if rs.get('warnings') or ts.get('warnings') else 'completed',
                         analysis_scope='recall_only_unverified_discovery_candidates',seconds=round(time.monotonic()-start,2),genome_bases=shards.total_bases,
                         candidate_pool=pool,validation=pool_validation,recall=rs,translated_search=ts,output=str(out))
            save(out/'sample_summary.json',summary);LOG.info('Recall-only completed sample=%s candidates=%d status=%s',args.sample,pool['candidate_model_records'],summary['status']);return summary
        if getattr(args,'annotation_mode','unified')=='unified':
            from .refine import annotate
            result=annotate(args,expanded,refs,evidence,out)
            reports,discovery_warnings=discovery_reports(out)
            if discovery_warnings:LOG.warning('[%s] Discovery incomplete: %s (records retained)',args.sample,', '.join(discovery_warnings))
            summary=dict(sample=args.sample,version=__version__,status='completed_with_warnings' if discovery_warnings or result['status']=='completed_with_warnings' else result['status'],seconds=round(time.monotonic()-start,2),genome_bases=shards.total_bases,
                analysis_scope='genome_discovery_structure_identity_annotation',candidate_pool=pool,candidate_validation=pool_validation,annotation=result,
                discovery=reports,discovery_warnings=discovery_warnings,output=str(out/'results'))
            save(out/'sample_summary.json',summary);return summary
        for m in union:
            cid='ORS2_'+args.sample+'_'+digest([m['contig'],m['strand'],m['exons'],m['cds'],m['frameshifts']])[:12]
            m.update(cid=cid,id=cid+'.original',source='original')
        evidence.apply(union,'union_'+args.sample)
        review_rows=[];review_models=[]
        for m in union:
            routing=route(m);review_rows.append(dict(candidate_id=m['cid'],contig=m['contig'],strand=m['strand'],start=min(s for s,e in m['exons']),end=max(e for s,e in m['exons']),
                length=m['length'],or_ga=int(m['or_ga']),or_evalue=m['or_evalue'],or_bits=m['or_bits'],raw_or_bits=m.get('raw_or_bits',0),gr_bits=m['gr_bits'],family=m['family'],family_support=m.get('family_support',''),route=routing,search_passes=','.join(m['search_passes'])))
            if routing=='OR_review':review_models.append(m)
        review.write_table(out/'discovery/candidate_pool.tsv',review_rows,list(review_rows[0]) if review_rows else ['candidate_id'])
        save(out/'discovery/summary.json',dict(first=len(first),second=len(second),translated=len(translated),union=len(union),first_only=sum(m['search_passes']==['first'] for m in union),
           routing=dict(Counter(r['route'] for r in review_rows)),adaptive_references=len(seeds),annotation_references=len(annotation),subga_review=sum(not m['or_ga'] for m in review_models)))
        LOG.info('OR structural review pool=%d (sub-GA=%d)',len(review_models),sum(not m['or_ga'] for m in review_models))
        rows,byid,decisions=structure(args,review_models,out,evidence,refs)
        tree_status='completed'
        try:flags=phylogeny(args,rows,decisions,refs,out)
        except Exception as error:
            LOG.exception('Phylogeny diagnostic failed; structural results will be exported with explicit missing-tree status')
            flags={};tree_status='failed'
            save(out/'phylogeny/diagnostic_failure.json',dict(error=str(error),structural_results_retained=True))
            for d in decisions.values():
                d['reason']+=';phylogeny_not_completed'
                if d['category']=='complete':d['confidence']='provisional_complete'
        # One bounded tree-triggered rescue for models that initially skipped local mapping.
        flagged=[r for r in rows if r['candidate_id'] in flags and (flags[r['candidate_id']]['tree_side']=='GR_side' or float(flags[r['candidate_id']]['gap_fraction'])>.5)]
        retry_records=[]
        for row in flagged:
            cid=row['candidate_id'];case=out/'structure/cases'/cid;status=json.loads((case/'local_status.json').read_text()) if (case/'local_status.json').exists() else {}
            if status.get('local_rebuilt'):retry_records.append(dict(candidate_id=cid,status='already_locally_rebuilt_unresolved'));continue
            original=next(m for m in byid[cid] if m['source']=='original');rids=choose_refs(original,refs,metadata=args.reference_metadata,policy=getattr(args,'reconstruction_reference_policy','audited'))
            retry=out/'structure/tree_retry';retry.mkdir(exist_ok=True)
            lengths=review.read_fai(Path(str(args.genome)+'.fai'))
            try:
                extra=initial_case(row,original,rids,retry,args,lengths,True);evidence.apply(extra,'tree_retry_'+cid)
                for m in extra:
                    if m['source']=='original':continue
                    m['fusion_pairs']=[];m['annotation_gene_ids']=[];byid[cid].append(m)
                # Do not elevate an un-crosschecked repair after a phylogenetic flag.
                retry_records.append(dict(candidate_id=cid,status='extra_models_saved_confidence_not_elevated'))
            except Exception as error:retry_records.append(dict(candidate_id=cid,status='error',error=str(error)))
        review.write_table(out/'structure/tree_retry.tsv',retry_records,list(retry_records[0]) if retry_records else ['candidate_id','status'])
        save_groups(out/'structure/all_models.json',byid)
        stats=dict(program='ORsearcher2',version=__version__,candidate_reviewed=len(review_models),phylogeny_flags=len(flagged),
          limitations=['Same assembly computational support, not experimental verification','Low-priority candidates are retained but not exhaustively structurally rebuilt','Reference-tree flags are not single-condition deletion rules'])
        if rows:
            export_args=SimpleNamespace(outdir=out)
            review.export_results(rows,byid,decisions,export_args,stats,{})
            unique=export_unique(rows,decisions,out,allow_containment=True)
            v1=validate(args.genome,out/'final');v2=validate(args.genome,out/'nonredundant')
        else:
            unique=dict(representative_loci=0,category_counts={});v1=v2=dict(checked_models=0)
            for folder in [out/'final',out/'nonredundant']:
                folder.mkdir(exist_ok=True);(folder/'OR_models.gff3').write_text('##gff-version 3\n');review.write_table(folder/'OR_catalogue.tsv',[],['candidate_id'])
                for k in ['complete','partial','pseudogene']:
                    write_fasta(folder/f'OR_{k}.proteins.faa',[]);write_fasta(folder/f'OR_{k}.cds.fna',[])
                save(folder/'OR_summary.json',unique)
                save(folder/'export_validation.json',dict(checked_models=0,coordinate_sequence_mismatches=0))
        cross_status=json.loads((out/'structure/crosscheck_summary.json').read_text());statuses=cross_status['counts']
        cross_warnings=[]
        if any(s['status']!='completed' for s in cross_status['statuses']):cross_warnings.append('secondary_alignment_unresolved_or_failed')
        if any(s.get('legacy_status')=='deferred_budget' for s in cross_status['statuses']):cross_warnings.append('legacy_reconstruction_budget_exceeded')
        translated_status=json.loads((out/'discovery/translated/summary.json').read_text())
        recall_status=json.loads((out/'discovery/recall/summary.json').read_text())
        summary=dict(sample=args.sample,status='completed_with_warnings' if tree_status!='completed' or translated_status.get('warnings') or recall_status.get('warnings') or cross_warnings else 'completed',phylogeny_status=tree_status,translated_search=translated_status,recall_search=recall_status,candidate_pool=pool,candidate_validation=pool_validation,crosscheck_warnings=cross_warnings,version=__version__,seconds=round(time.monotonic()-start,2),genome_bases=shards.total_bases,search_records=len(union),
            review_pool=len(review_models),family_unresolved=len(review_models)-len(rows),representatives=unique,exonerate_statuses=statuses,validation=dict(all_records=v1,representatives=v2),output=str(out))
        save(out/'sample_summary.json',summary);LOG.info('Sample completed %s',summary);return summary
    finally:LOG.removeHandler(handler);handler.close()

def run(args):
    root=args.outdir.resolve()
    if root.exists() and any(root.iterdir()) and not args.resume:raise FileExistsError('New run directory required; use --resume only for identical inputs/config')
    root.mkdir(parents=True,exist_ok=True)
    from .logging_utils import configure, log_start
    configure(root,getattr(args,'log_level','info'),getattr(args,'quiet',False))
    # Native tools with an explicit thread option receive the requested
    # budget. Hidden numerical/OpenMP pools must not multiply local workers.
    for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
        os.environ[variable]='1'
    samples=load_samples(args)
    if not samples or len({s['sample'] for s in samples})!=len(samples):raise ValueError('No samples or duplicate sample IDs')
    log_start(LOG, __version__, args, root, len(samples))
    additional=([args.search_seeds] if getattr(args,'search_seeds',None) else [])+([args.complete_reference] if getattr(args,'complete_reference',None) else [])
    if getattr(args,'previous_sample',None):
        additional.extend(args.previous_sample/name for name in ['structure/all_models.json','annotation_qc/identity_audit.tsv','nonredundant/OR_catalogue.tsv','structure/crosscheck_summary.json'] if (args.previous_sample/name).is_file())
    inputs=[args.or_reference,args.gr_reference,args.or_hmm,*args.gr_hmm]+([args.exclude_reference_proteins] if args.exclude_reference_proteins else [])+[Path(s['genome']) for s in samples]+[Path(s['annotation']) for s in samples if s.get('annotation')]+getattr(args,'candidate_models',[])+([args.background_reference] if getattr(args,'background_reference',None) else [])+additional
    source_hash=digest([(str(p.relative_to(Path(__file__).parent)),file_digest(p)) for p in sorted(Path(__file__).parent.rglob('*')) if p.suffix in ['.py','.R']])
    small_reference_inputs=[args.or_reference,args.gr_reference,args.or_hmm,*args.gr_hmm]+([args.background_reference] if getattr(args,'background_reference',None) else [])+([args.complete_reference] if getattr(args,'complete_reference',None) else [])
    extra_inputs = getattr(args,'extra_reference_inputs',[])
    fingerprint=dict(version=__version__,source_sha256=source_hash,reference_sha256={str(p):file_digest(p) for p in [*small_reference_inputs,*extra_inputs]},inputs=[metadata(p) for p in [*inputs,*extra_inputs]],parameters={k:v for k,v in vars(args).items() if k not in ['resume','func','log_level','quiet']})
    fingerprint=json.loads(json.dumps(fingerprint,default=str))
    if (root/'input_manifest.json').exists() and json.loads((root/'input_manifest.json').read_text())!=fingerprint:raise ValueError('Resume inputs/config changed')
    save(root/'input_manifest.json',fingerprint)
    from .tuning import policy
    save(root/'effective_parameters.json',dict(parameters={k:v for k,v in vars(args).items() if k!='thresholds'},thresholds=policy(args=args)))
    LOG.info('Checking analysis tools and preparing reference inputs')
    versions=dependency_check(args);save(root/'dependencies.json',dict(cpu_logical=os.cpu_count(),tools=versions))
    from .input_files import materialize_references
    materialize_references(args, root)
    refs=prepare_references(args,root);evidence=Evidence(args,root/'shared/evidence',refs);results=[]
    if getattr(args,'evidence_cache',None):
        import sqlite3
        cache=Path(args.evidence_cache);summary=json.loads((cache.parent/'summary.json').read_text())
        if summary['config']!=evidence.config:raise ValueError('Incompatible evidence cache; references/HMM/audit differ')
        old=sqlite3.connect(cache.resolve().as_uri()+'?mode=ro',uri=True);old.backup(evidence.connection);old.close()
        save(root/'shared/evidence/cache_source.json',dict(path=str(cache),read_only_source=True,config=evidence.config))
    for index, spec in enumerate(samples, 1):
        local=SimpleNamespace(**vars(args));local.sample=spec['sample'];local.genome=Path(spec['genome']);local.annotation=Path(spec['annotation']) if spec.get('annotation') else None
        LOG.info('Sample %d/%d | %s', index, len(samples), local.sample)
        try:results.append(run_sample(local,root,refs,evidence))
        except Exception as error:
            LOG.exception('Sample failed %s',local.sample);record=dict(sample=local.sample,status='error',error=str(error));results.append(record);save(root/'samples'/local.sample/'sample_summary.json',record)
        save(root/'batch_summary.json',dict(version=__version__,samples=results,status='running'))
    failed=any(r['status']=='error' for r in results)
    warnings=any(r['status']=='completed_with_warnings' for r in results)
    status='failed' if failed else 'completed_with_warnings' if warnings else 'completed'
    save(root/'batch_summary.json',dict(version=__version__,samples=results,status=status))
    evidence.connection.close()
    LOG.info('Batch %s | samples=%d | output=%s',status,len(results),root)
    return 1 if failed else 0
