"""Lossless shared sequence/evidence payloads for structural model groups.

Genomic placements and true aligner provenance stay distinct. This is storage
deduplication, never genomic locus merging. Legacy dict-of-lists remains readable.
"""
import copy,json
from pathlib import Path
from .runtime import digest,save

PAYLOAD_FIELDS={'protein','cds','homology_hits','raw_homology_hits','gr_homology_hits',
    'or_profile_domains','gr_profile_domains','reference_exact_matches',
    'qualified_reference_exact_matches','full_ref_support','background_hits'}


def pack_groups(groups):
    payloads={};packed={}
    for cid,models in groups.items():
        values=[]
        for m in models:
            payload={k:v for k,v in m.items() if k in PAYLOAD_FIELDS}
            key=digest(payload);payloads.setdefault(key,payload)
            placement={k:v for k,v in m.items() if k not in PAYLOAD_FIELDS and k!='discovery_alternatives'}
            if m.get('discovery_alternatives'):
                placement['discovery_alternative_ids']=[a['id'] for a in m['discovery_alternatives']]
            placement['_sequence_evidence_ref']=key;values.append(placement)
        packed[cid]=values
    return dict(format='orsearcher-model-groups-v1',payloads=payloads,groups=packed)


def save_groups(path,groups):
    save(path,pack_groups(groups),compact=True)


def load_groups(path):
    data=json.loads(Path(path).read_text(encoding='utf-8'))
    if data.get('format')!='orsearcher-model-groups-v1':return data
    groups={}
    for cid,placements in data['groups'].items():
        groups[cid]=[]
        for p in placements:
            m=copy.deepcopy({k:v for k,v in p.items() if k!='_sequence_evidence_ref'})
            m.update(copy.deepcopy(data['payloads'][p['_sequence_evidence_ref']]))
            groups[cid].append(m)
    return groups
