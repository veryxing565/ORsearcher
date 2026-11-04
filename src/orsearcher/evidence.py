from __future__ import annotations
from .tuning import policy as tuning_policy, reconstruction_metadata
import json,math,sqlite3
from collections import defaultdict
from pathlib import Path
from .fasta import FastaRecord,write_fasta
from .runtime import digest,execute,file_digest,save,LOG

def ga_thresholds(hmm):
    from .input_files import read_profiles
    return {model.name: model.ga for model in read_profiles(hmm)}


def profile_cutoff(domain, cutoffs):
    return cutoffs.get(domain.get('profile')) if isinstance(cutoffs, dict) else cutoffs


def passes_ga(domain, cutoffs):
    cutoff = profile_cutoff(domain, cutoffs)
    return cutoff is not None and domain['sequence_score'] >= cutoff[0] and domain['domain_score'] >= cutoff[1]

def parse_domains(path):
    values=defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split();values[f[0]].append(dict(profile=f[3],sequence_evalue=float(f[6]),sequence_score=float(f[7]),domain_evalue=float(f[12]),domain_score=float(f[13]),
            coverage=(int(f[16])-int(f[15])+1)/int(f[5]),start=int(f[17]),end=int(f[18])))
    return values

def alignment_provenance(model,metadata):
    """Self/annotation queries aid discovery but never independent validation."""
    ref=model.get('ref','');info=metadata.get(ref,{})
    origin=('self_generated' if ref.startswith('adaptive_') else 'target_annotation' if ref.startswith('annotation_')
            else 'external' if ref in metadata else 'unknown')
    model.update(alignment_reference_origin=origin,
                 alignment_reference_qualified=origin=='external' and bool(info.get('reconstruction_allowed',False)),
                 reference_sequence_key=info.get('sequence_key',''))
    return model

class Evidence:
    """Protein evidence shared by exact sequence; never merges genomic loci."""
    def __init__(self,args,folder,refs):
        self.args=args;self.folder=Path(folder);self.folder.mkdir(parents=True,exist_ok=True);self.refs=refs
        self.connection=sqlite3.connect(self.folder/'proteins.sqlite')
        self.connection.execute('CREATE TABLE IF NOT EXISTS evidence (key TEXT PRIMARY KEY, data TEXT NOT NULL)')
        self.reference_metadata=getattr(args,'reference_metadata',{})
        self.sequence_references=defaultdict(list)
        for rid,r in refs.items():self.sequence_references[r.sequence].append(rid)
        self.config=digest(dict(version='2.1.0-evidence-2',reference=file_digest(args.or_reference),gr_reference=file_digest(Path(str(args.gr_db)+'.dmnd')),reference_audit=self.reference_metadata,
             hmms=[file_digest(p) for p in [args.or_hmm,*args.gr_hmm]],commands='diamond sensitive top50OR top5GR, hmm domtbl E10',tuning=tuning_policy(args=args)))
        self.hmm_ga=[ga_thresholds(p) for p in [args.or_hmm,*args.gr_hmm]];self.batch=0;self.computed=0;self.reused=0

    def apply(self,models,label):
        pending={};entries={}
        for m in models:
            protein=m['protein'].replace('*','X');key=digest([self.config,protein]);m['evidence_key']=key
            if key in entries:continue
            r=self.connection.execute('SELECT data FROM evidence WHERE key=?',(key,)).fetchone()
            if r:entries[key]=json.loads(r[0]);self.reused+=1
            elif protein:pending[key]=protein
            else:entries[key]=dict(hmms=[[] for _ in self.hmm_ga],or_hits=[],gr_hits=[])
        if pending:
            self.batch+=1;batch=self.folder/f'{self.batch:04d}_{label}';batch.mkdir(exist_ok=True)
            write_fasta(batch/'proteins.faa',[FastaRecord(k,'',p) for k,p in pending.items()])
            hmmhits=[]
            for i,hmm in enumerate([self.args.or_hmm,*self.args.gr_hmm]):
                dom=batch/f'hmm{i}.domtbl'
                execute(['hmmsearch','--cpu',self.args.threads,'--noali','-E','10','--domE','10','--domtblout',dom,hmm,batch/'proteins.faa'],batch/f'hmm{i}.out',inputs=[hmm,batch/'proteins.faa'])
                hmmhits.append(parse_domains(dom))
            families={}
            for family,db,limit in [('or',self.args.or_db,getattr(self.args,'evidence_or_targets',50)),('gr',self.args.gr_db,getattr(self.args,'evidence_gr_targets',5))]:
                table=batch/f'{family}.tsv'
                execute(['diamond','blastp','--query',batch/'proteins.faa','--db',db,'--out',table,'--outfmt','6','qseqid','sseqid','pident','length','qlen','slen','qstart','qend','sstart','send','evalue','bitscore',
                    '--max-target-seqs',limit,'--sensitive','--evalue',getattr(self.args,'evidence_evalue',1e-3),'--threads',self.args.threads,'--quiet'],batch/f'{family}.command.log',inputs=[batch/'proteins.faa',Path(str(db)+'.dmnd')])
                hits=defaultdict(list)
                for line in table.read_text().splitlines():
                    f=line.split('\t');hits[f[0]].append(dict(ref=f[1],identity=float(f[2]),length=int(f[3]),qlen=int(f[4]),slen=int(f[5]),qs=int(f[6]),qe=int(f[7]),ss=int(f[8]),se=int(f[9]),evalue=float(f[10]),bits=float(f[11])))
                families[family]=hits
            for key in pending:
                entry=dict(hmms=[h.get(key,[]) for h in hmmhits],or_hits=families['or'].get(key,[]),gr_hits=families['gr'].get(key,[]))
                entries[key]=entry;self.connection.execute('INSERT OR REPLACE INTO evidence VALUES (?,?)',(key,json.dumps(entry)))
            self.connection.commit();self.computed+=len(pending)
            LOG.debug('Evidence %s: models=%d new_unique=%d',label,len(models),len(pending))
        for m in models:self.assign(m,entries[m['evidence_key']])
        save(self.folder/'summary.json',dict(new_unique_proteins=self.computed,cache_hits=self.reused,config=self.config))

    def assign(self,m,e):
        from .tuning import policy
        t=policy(args=getattr(self,'args',None));m['_thresholds']=dict(t);m['_evidence_config']=getattr(self,'config',None)
        domains=e['hmms'][0];oh=e['or_hits'];gh=e['gr_hits'];best=min(domains,key=lambda x:x['sequence_evalue']) if domains else {}
        og=any(passes_ga(d,self.hmm_ga[0]) for d in domains)
        gg=any(passes_ga(d,ga) for ds,ga in zip(e['hmms'][1:],self.hmm_ga[1:]) for d in ds)
        from .references import supported_overlap
        qualified=[]
        for hit in oh:
            info=self.reference_metadata.get(hit['ref'],{})
            hit['domain_overlap']=supported_overlap(hit,info.get('domains',[]))
            hit['curated']=info.get('curated',False)
            hit['reference_sequence_key']=info.get('sequence_key',hit['ref'])
            # A profile-negative curated reference remains usable, but an
            # automatic reference's non-OR extension is not family evidence.
            if hit['domain_overlap']>=t['evidence_domain_overlap'] or info.get('curated',False):qualified.append(hit)
        ob=qualified[0]['bits'] if qualified else 0.;gb=gh[0]['bits'] if gh else 0.
        curated_keys={h['reference_sequence_key'] for h in qualified if h['curated'] and h['bits']>=100 and h['evalue']<=1e-5
                      and (h['qe']-h['qs']+1)/h['qlen']>=.6 and h['identity']>=25}
        domain_signal=og or best.get('sequence_evalue',1.)<=1e-3
        family_support=domain_signal or len(curated_keys)>=2
        m.update(or_evalue=best.get('sequence_evalue',1.),or_score=best.get('sequence_score',0.),or_ga=og,gr_ga=gg,or_bits=ob,gr_bits=gb,
                 or_profile_domains=domains,gr_profile_domains=e['hmms'][1:],gr_homology_hits=gh,
                 family='GR_likely' if gb-ob>=10 and gb else 'OR' if qualified and ob-gb>=10 and family_support else 'OR_GR_uncertain',homology_hits=qualified,
                 raw_or_bits=oh[0]['bits'] if oh else 0.,raw_homology_hits=oh,curated_reference_count=len(curated_keys),
                 family_support='OR_profile' if domain_signal else 'curated_homology_consensus' if len(curated_keys)>=2 else 'insufficient',
                 full_ref_support=[],reference_exact_matches=list(self.sequence_references.get(m['protein'],[])),repeated_or_domains=False,repeated_reference_regions=False)
        alignment_provenance(m,self.reference_metadata)
        m['qualified_reference_exact_matches']=[r for r in m['reference_exact_matches'] if self.reference_metadata.get(r,{}).get('reconstruction_allowed',False)]
        full_sequence_keys=set()
        for hit in qualified:
            ref=self.refs.get(hit['ref'])
            slack_s=max(20,math.ceil(hit['slen']*.1));slack_q=max(20,math.ceil(m['length']*.1))
            if (ref and self.reference_metadata.get(hit['ref'],{}).get('reconstruction_allowed',False) and hit['evalue']<=1e-5 and hit['identity']>=25 and (hit['se']-hit['ss']+1)/hit['slen']>=.85
                and hit['ss']<=slack_s and hit['slen']-hit['se']<=slack_s and .75<=m['length']/hit['slen']<=1.25
                and hit['qs']<=slack_q and m['length']-hit['qe']<=slack_q):
                sequence_key=hit['reference_sequence_key']
                if sequence_key not in full_sequence_keys:
                    m['full_ref_support'].append(hit['ref']);full_sequence_keys.add(sequence_key)
        eligible=[d for d in domains if d['domain_evalue']<=1e-3 and d['coverage']>=.5
                  and d['domain_score']>=(profile_cutoff(d,self.hmm_ga[0]) or (0.,0.))[1]]
        m['repeated_or_domains']=any(a.get('profile')==b.get('profile') and min(a['end'],b['end'])-max(a['start'],b['start'])+1<=5
                                    for i,a in enumerate(eligible) for b in eligible[i+1:])
        if getattr(getattr(self,'args',None),'identity_policy','legacy')=='multichannel':
            from .identity import structural_family
            structural_family(m)

def route(model):
    """No hard OR-GA gate: sub-GA and strong homologs enter structural review."""
    strong=model.get('curated_reference_count',0)>=2 and model['or_bits']>=100
    subga=model['or_evalue']<=1e-3 and model['or_bits']>=50
    if model['family']=='GR_likely' and model['gr_ga'] and not model['or_ga'] and not subga:return 'GR_likely'
    if model['or_ga'] or subga or strong:return 'OR_review'
    return 'low_priority'
