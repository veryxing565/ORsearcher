"""Audit full references; mask unsupported regions of automatic gap-fill queries.

Original references remain intact for provenance and local reconstruction.
Masking changes discovery queries, never the genome or an exported CDS.
"""
from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from .fasta import FastaRecord, read_fasta, write_fasta
from .evidence import parse_domains, ga_thresholds, passes_ga, profile_cutoff
from .runtime import execute, save, digest


def supported_overlap(hit, domains):
    intervals = sorted((max(hit['ss'], d['start']), min(hit['se'], d['end']))
                       for d in domains if d['domain_evalue'] <= 1e-3)
    merged = []
    for left, right in intervals:
        if right < left:
            continue
        if merged and left <= merged[-1][1] + 1:
            merged[-1][1] = max(right, merged[-1][1])
        else:
            merged.append([left, right])
    return sum(b-a+1 for a, b in merged) / max(1, hit['se']-hit['ss']+1)


def audit_references(args, folder, records):
    folder = Path(folder)
    path = folder/'or.faa'
    dom = folder/'or_reference.domtbl'
    execute(['hmmsearch', '--cpu', args.threads, '--noali', '-E', '10',
             '--domE', '10', '--domtblout', dom, args.or_hmm, path],
            folder/'or_reference_hmm.log', inputs=[args.or_hmm, path])
    domains = parse_domains(dom)
    meta = {}
    queries = []
    ga = ga_thresholds(args.or_hmm)
    user_ids = getattr(args,'user_or_ids',set())
    for r in records:
        ds = [d for d in domains.get(r.id, []) if d['domain_evalue'] <= 1e-3]
        automatic = r.id not in user_ids and bool(re.search(r'uniprotgap|gap_fill|self_generated', r.description, re.I))
        curated = r.id not in user_ids and bool(re.search(r'source=(?:iOrME|Mitchell)|^iorme\||^mitchell\|', r.description, re.I))
        profile_counts = Counter(d['profile'] for d in ds)
        flags = []
        if not ds:
            flags.append('no_OR_domain')
        if any(count>1 for count in profile_counts.values()):
            flags.append('multiple_OR_domain_segments')
        if ds and max(d['coverage'] for d in ds) < .75:
            flags.append('incomplete_profile_coverage')
        if not 250 <= len(r.sequence) <= 650:
            flags.append('unusual_length')
        sequence = r.sequence
        # Automatic taxonomic-gap references must not seed thousands of loci
        # solely through a domain-negative extension. Keep profile coordinates.
        if automatic:
            keep = [False] * len(sequence)
            for d in ds:
                for pos in range(max(0, d['start']-1), min(len(sequence), d['end'])):
                    keep[pos] = True
            sequence = ''.join(aa if keep[i] else 'X' for i, aa in enumerate(sequence))
        reconstruction = (250 <= len(r.sequence) <= 650 and (curated or any(
            profile_counts[d['profile']] == 1 and d['coverage'] >= .75 and
            (profile_cutoff(d,ga) is None or passes_ga(d,ga)) for d in ds)))
        meta[r.id] = dict(curated=curated, automatic=automatic, domains=ds,
                         reconstruction_allowed=reconstruction,
                         sequence_key=digest(r.sequence), flags=flags)
        if not automatic or any(aa != 'X' for aa in sequence):
            queries.append(FastaRecord(r.id, r.description, sequence))
    write_fasta(folder/'or.discovery.faa', queries)
    save(folder/'reference_audit.json', meta)
    save(folder/'reference_audit_summary.json', dict(input_references=len(records),
         discovery_references=len(queries), automatic_masked=sum(m['automatic'] for m in meta.values()),
         reconstruction_references=sum(m['reconstruction_allowed'] for m in meta.values()),
         policy='Full originals retained; automatic gap-fill discovery queries masked outside OR domain; reconstruction separately qualified'))
    args.reference_metadata = meta
    args.discovery_references = {r.id: r for r in queries}
    return meta
