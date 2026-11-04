#!/usr/bin/env python3
"""Build the pinned public panels from archived UniProtKB search TSV files."""
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile

from orsearcher.evidence import ga_thresholds, parse_domains, passes_ga
from orsearcher.fasta import FastaRecord, write_fasta
from orsearcher.input_files import PROTEIN


def candidates(snapshot, sources, family):
    rows = {}
    for source in sources:
        if source['family'] != family:
            continue
        with (snapshot/source['file']).open(encoding='utf-8',newline='') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                seq = row['Sequence'].rstrip('*')
                if 100 <= len(seq) <= 700 and not (set(seq)-PROTEIN):
                    rows[row['Entry']] = (row,source['order'],seq)
    return rows


def screen(rows, models, work, threads):
    target = work/'proteins.faa'
    write_fasta(target, (FastaRecord(rid,'',entry[2]) for rid,entry in sorted(rows.items())))
    all_hits = defaultdict(list)
    cutoffs = {}
    for index, model in enumerate(models):
        output = work/f'model_{index}.domtbl'
        subprocess.run(['hmmsearch','--cpu',str(threads),'--noali','--domE','0.001',
                        '--domtblout',str(output),str(model),str(target)],
                       check=True,stdout=subprocess.DEVNULL)
        cutoffs.update(ga_thresholds(model))
        for rid,hits in parse_domains(output).items():
            all_hits[rid].extend(d for d in hits if d['domain_evalue']<=.001 and d['coverage']>=.25)
    return all_hits,cutoffs


def select(rows, hits, cutoffs, family, exclude):
    by_order,by_taxon = Counter(),Counter()
    seen = set(exclude)
    kept,complete = [],[]
    keys = sorted((rid for rid in rows if hits.get(rid)), key=lambda rid:(
        rows[rid][0]['Reviewed']!='reviewed', min(d['domain_evalue'] for d in hits[rid]), rid))
    for rid in keys:
        row,order,seq = rows[rid]
        taxon = row['Organism (ID)']
        if seq in seen or by_taxon[taxon]>=96 or by_order[order]>=250:
            continue
        rec = FastaRecord(f'uniprot_{family}_{rid}',f'source=UniProtKB accession={rid} taxid={taxon} order={order}',seq)
        kept.append(rec)
        if (family=='OR' and 250<=len(seq)<=650 and 'fragment' not in row['Protein names'].lower()
                and len(hits[rid])==1 and hits[rid][0]['coverage']>=.75 and passes_ga(hits[rid][0],cutoffs)):
            complete.append(rec)
        seen.add(seq)
        by_taxon[taxon]+=1
        by_order[order]+=1
    return kept,complete,dict(by_order)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,required=True,help='Directory with source.json and the referenced UniProt TSVs')
    parser.add_argument('--data',type=Path,required=True,help='Release bundle directory containing the pinned Pfam profiles and background')
    parser.add_argument('--threads',type=int,default=4)
    args = parser.parse_args()
    data=args.data
    sources=json.loads((args.snapshot/'source.json').read_text())
    if any('error' in source for source in sources):
        raise ValueError('At least one UniProt snapshot download failed')
    files=json.loads((data/'manifest.json').read_text())['files']
    license_info=dict(license='CC-BY-4.0',attribution='UniProt Consortium',
                      terms='https://www.uniprot.org/help/license',
                      license_url='https://creativecommons.org/licenses/by/4.0/')
    retained={}
    for family, names in (('OR',['OR.hmm']),('GR',['GR_1.hmm','GR_2.hmm'])):
        rows=candidates(args.snapshot,sources,family)
        with tempfile.TemporaryDirectory() as temp:
            hits,cutoffs=screen(rows,[data/n for n in names],Path(temp),args.threads)
        kept,complete,orders=select(rows,hits,cutoffs,family,{r.sequence for r in retained.get('OR',[])})
        if len(kept)<100 or len(orders)<5 or (family=='OR' and not complete):
            raise ValueError(f'Unexpectedly sparse UniProt {family} panel: {orders}')
        retained[family]=kept
        name=family+'.faa'
        write_fasta(data/name,kept)
        files[name]=dict(**license_info,records=len(kept),orders=orders,role=family+' protein reference',
            modification='Name and HMM selected, amino-acid filtered, sequence deduplicated, capped per taxon and order, and reheadered.')
        if family=='OR':
            write_fasta(data/'complete_OR.faa',complete)
            files['complete_OR.faa']=dict(**license_info,records=len(complete),role='Computationally qualified OR template subset',
                modification='Subset of OR.faa with 250-650 aa, one supported OR domain, at least 75% profile coverage, Pfam GA and no fragment label.')
    for name,item in files.items():
        blob=(data/name).read_bytes()
        item.update(bytes=len(blob),sha256=sha256(blob).hexdigest())
    for source in sources:
        source['sha256']=sha256((args.snapshot/source['file']).read_bytes()).hexdigest()
    manifest=dict(schema='orsearcher-reference-panel-v2',panel_id='public-uniprotkb-2026_03-pfam-2.3.0',
        files=files,uniprot_queries=sources,
        selection=dict(domain_evalue_max=.001,profile_coverage_min=.25,per_taxon_cap=96,per_order_cap=250),
        panel_scope='Named UniProtKB entries from eight sampled orders. See file-level coverage. The panel is not a held-out benchmark.')
    (data/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({name:{k:v for k,v in item.items() if k in ('records','orders','sha256')} for name,item in files.items()},indent=2))


if __name__=='__main__':
    main()
