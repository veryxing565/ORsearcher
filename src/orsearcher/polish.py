"""Practical OR locus consolidation and bounded, assembly-only repair.

This is a reproducible computational selection policy, not a calibrated
probability of correctness. Alternatives, failures and unresolved identities
are retained in the audit; no bases are grafted from a reference protein.
"""
from __future__ import annotations
from .tuning import policy as tuning_policy, reconstruction_metadata

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import re
import subprocess
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

from .fasta import FastaRecord, read_fasta, write_fasta
from .identity import background_evidence, judge
from .runtime import LOG, digest, metadata, save


def table(path):
    path = Path(path)
    if not path.exists():
        return []
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_table(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else ['model_id', 'reason'])
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def coding_length(model):
    return sum(e - s + 1 for s, e in model['exons'])


def bounds(model):
    return min(s for s, e in model['exons']), max(e for s, e in model['exons'])


def placement(model):
    return model['contig'], model['strand'], tuple(sorted(map(tuple, model['exons']))), model['cds']


def overlap(a, b):
    shared = sum(max(0, min(e, v) - max(s, u) + 1)
                 for s, e in a['exons'] for u, v in b['exons'])
    return shared / max(1, coding_length(a)), shared / max(1, coding_length(b)), shared


def junctions(model):
    exons = sorted(model['exons'])
    return {(a[1] + 1, b[0] - 1) for a, b in zip(exons, exons[1:])}


def same_frame(a, b):
    if a['strand'] != b['strand']:
        return False
    def offsets(model):
        n = 0
        result = []
        for start, end in sorted(model['exons'], reverse=model['strand'] == '-'):
            result.append((start, end, n))
            n += end - start + 1
        return result
    found = False
    for s, e, n in offsets(a):
        for t, u, k in offsets(b):
            if max(s, t) > min(e, u):
                continue
            found = True
            pos = max(s, t) if a['strand'] == '+' else min(e, u)
            x = n + (pos - s if a['strand'] == '+' else e - pos)
            y = k + (pos - t if b['strand'] == '+' else u - pos)
            if (x - y) % 3:
                return False
    return found


def splice_compatible(a, b):
    for left, right in ((a, b), (b, a)):
        exons = sorted(left['exons'])
        for x, y in zip(exons, exons[1:]):
            if any(min(y[0] - 1, e) >= max(x[1] + 1, s) for s, e in right['exons']):
                return False
    return True


def same_locus(a, b):
    """Fixed coding overlap/structure, never span proximity or peptide alone."""
    if a['contig'] != b['contig'] or a['strand'] != b['strand']:
        return False
    x, y, shared = overlap(a, b)
    if shared < 90:
        return False
    common_splice = bool(junctions(a) & junctions(b))
    if min(x, y) >= .65 and common_splice:
        return True
    if min(x, y) >= .85 and same_frame(a, b) and splice_compatible(a, b):
        return True
    small_cov, big_cov = (x, y) if coding_length(a) <= coding_length(b) else (y, x)
    return bool(small_cov >= .85 and big_cov >= .12 and
                (common_splice or (same_frame(a, b) and splice_compatible(a, b)
                 and (a['protein'] in b['protein'] or b['protein'] in a['protein']))))


def complete_reference_ids(args, refs):
    """Explicit complete-reference set overrides audited reconstruction IDs."""
    if hasattr(args, 'complete_reference_ids'):
        return set(args.complete_reference_ids) & set(refs), 'active_audited_complete_reference_ids_after_holdout'
    explicit = getattr(args, 'complete_reference', None)
    if explicit:
        records = list(read_fasta(Path(explicit)))
        missing = [r.id for r in records if r.id not in refs or refs[r.id].sequence != r.sequence]
        if missing:
            raise ValueError('Complete references must match the active OR reference IDs and sequences: ' + ', '.join(missing[:5]))
        return {r.id for r in records}, 'explicit_complete_reference'
    metadata = getattr(args, 'reference_metadata', {})
    return {rid for rid in refs if metadata.get(rid, {}).get('reconstruction_allowed', False)}, 'audited_reconstruction_eligible_not_experimental_completeness'


def assessment(model, core_ids):
    from .tuning import policy
    t=policy(model)
    hits = [h for h in model.get('homology_hits', [])
            if h['evalue'] <= t['model_hit_evalue'] and h['bits'] >= t['model_hit_bits'] and h['identity'] >= t['model_hit_identity']
            and (h['qe'] - h['qs'] + 1) / max(1, h['qlen']) >= t['model_hit_coverage']]
    whole = [h for h in hits if (h['ref'] in core_ids or t['complete_evidence']=='external-homology') and not h['ref'].startswith(('adaptive_','annotation_')) and h['bits'] >= t['complete_hit_bits']
             and (h['qe'] - h['qs'] + 1) / max(1, h['qlen']) >= t['complete_query_coverage']
             and (h['se'] - h['ss'] + 1) / max(1, h['slen']) >= t['complete_reference_coverage']
             and h['ss'] <= max(t['complete_terminal_slack_aa'], t['complete_terminal_slack_fraction'] * h['slen'])
             and h['slen'] - h['se'] <= max(t['complete_terminal_slack_aa'], t['complete_terminal_slack_fraction'] * h['slen'])]
    best = max(whole or hits, key=lambda h: (h['bits'], h['identity']), default=None)
    identity = judge(model)
    identity_ok = identity['status'] == 'OR_supported'
    intact = not any(model.get(k) for k in ('frameshifts', 'stops', 'ambiguous', 'remainder'))
    safe = (bool(model.get('canonical')) and not model.get('bridge') and not model.get('fusion_pairs')
            and not (model.get('repeated_or_domains') and model.get('repeated_reference_regions')))
    complete = bool(identity_ok and safe and intact and model.get('has_start') and model.get('has_stop')
                    and t['complete_min_aa'] <= model['length'] <= t['complete_max_aa'] and (whole or t['complete_evidence']=='coding-and-identity'))
    partial = bool(identity_ok and safe and intact and t['partial_min_aa'] <= model['length'] <= t['partial_max_aa']
                   and any(h['identity'] >= t['partial_hit_identity'] for h in hits))
    blockers = [name for name, passed in [
        ('OR_identity_unresolved', identity_ok), ('unsafe_structure', safe),
        ('coding_damage_or_ambiguity', intact), ('missing_start', model.get('has_start')),
        ('missing_stop', model.get('has_stop')), ('length_outside_complete_range', t['complete_min_aa'] <= model['length'] <= t['complete_max_aa']),
        ('no_qualified_full_length_template_support', bool(whole) or t['complete_evidence']=='coding-and-identity')] if not passed]
    return dict(usable_complete=complete, credible_partial=partial, identity_ok=identity_ok,
                completeness_blockers=blockers,completeness_policy=t['complete_evidence'],full_length_template_supported=any(h['ref'] in core_ids for h in whole),
                identity_status=identity['status'], identity_reason=identity['reason'],
                intact=intact, safe=safe, best_reference=best['ref'] if best else '',
                core_reference=bool(best and best['ref'] in core_ids),
                reference_coverage=(best['se'] - best['ss'] + 1) / best['slen'] if best else 0,
                query_coverage=(best['qe'] - best['qs'] + 1) / best['qlen'] if best else 0,
                identity_pct=best['identity'] if best else 0, bits=best['bits'] if best else 0)


def rank(model):
    a = model['_assessment']
    return (a['usable_complete'], a['credible_partial'], a['identity_ok'],
            model.get('_reviewed_complete', False), a['core_reference'],
            round(a['reference_coverage'], 2), a['bits'], a['identity_pct'],
            model.get('_baseline_selected', False),
            bool(model.get('has_start') and model.get('has_stop')), model['id'])


def consolidate(models):
    """Non-transitive anchors; retain competing/bridging models in audit."""
    clusters, assignments, conflicts = [], [], []
    bins, strong_bins = defaultdict(set), defaultdict(set)
    ordered = sorted(models, key=rank, reverse=True)
    strong = []
    for m in ordered:
        a = m['_assessment']
        if not (a['usable_complete'] and a['reference_coverage'] >= .9 and a['query_coverage'] >= .9 and a['identity_pct'] >= 40):
            continue
        left, right = bounds(m)
        near = {i for b in range(left // 100000, right // 100000 + 1) for i in strong_bins[(m['contig'], m['strand'], b)]}
        if any(same_locus(m, strong[i]) for i in near):
            continue
        index = len(strong)
        strong.append(m)
        for b in range(left // 100000, right // 100000 + 1):
            strong_bins[(m['contig'], m['strand'], b)].add(index)
    for m in ordered:
        if not (m['_assessment']['usable_complete'] or m['_assessment']['credible_partial']):
            assignments.append(dict(model_id=m['id'], representative_model_id='', disposition='not_usable_main_set', reason=m['_assessment']['identity_reason']))
            continue
        left, right = bounds(m)
        nearby = {i for b in range(left // 100000, right // 100000 + 1) for i in strong_bins[(m['contig'], m['strand'], b)]}
        anchors = [strong[i] for i in nearby if overlap(m, strong[i])[1] >= .35 and overlap(m, strong[i])[2] >= 150]
        if any(bounds(a)[1] < bounds(b)[0] or bounds(b)[1] < bounds(a)[0] for i, a in enumerate(anchors) for b in anchors[i + 1:]):
            assignments.append(dict(model_id=m['id'], representative_model_id='', disposition='possible_two_gene_fusion', reason='bridges_separate_strong_coding_anchors'))
            continue
        nearby = {i for b in range(left // 100000, right // 100000 + 1) for i in bins[(m['contig'], m['strand'], b)]}
        matches = [i for i in sorted(nearby) if same_locus(m, clusters[i][0])]
        if len(matches) > 1:
            scores = sorted(((overlap(m, clusters[i][0])[2], i) for i in matches), reverse=True)
            if scores[0][0] < 1.5 * scores[1][0]:
                assignments.append(dict(model_id=m['id'], representative_model_id='', disposition='ambiguous_multi_anchor_overlap', reason='no_transitive_merge'))
                continue
            matches = [scores[0][1]]
        if matches:
            index = matches[0]
            clusters[index].append(m)
            disposition = 'same_locus_alternative'
        else:
            index = len(clusters)
            clusters.append([m])
            for b in range(left // 100000, right // 100000 + 1):
                bins[(m['contig'], m['strand'], b)].add(index)
            disposition = 'representative'
        assignments.append(dict(model_id=m['id'], representative_model_id=clusters[index][0]['id'], disposition=disposition, reason='fixed_coding_anchor'))
    # Practical one-representative policy is intentionally distinct from an
    # assertion that every overlapping path is the same biological gene.
    blocked = set()
    for i, cluster in enumerate(clusters):
        if i in blocked:
            continue
        a = cluster[0]
        for j in range(i + 1, len(clusters)):
            if j in blocked:
                continue
            b = clusters[j][0]
            if a['contig'] != b['contig'] or a['strand'] != b['strand']:
                continue
            x, y, shared = overlap(a, b)
            if shared >= 30:
                blocked.add(j)
                conflicts.append(dict(retained_model=a['id'], held_model=b['id'], reason='competing_coding_overlap_not_independent_main_gene', overlap_nt=shared, overlap_retained=x, overlap_held=y))
    return [c for i, c in enumerate(clusters) if i not in blocked], assignments, conflicts


def load_pool(source, groups, decisions, core_ids, thresholds=None):
    source = Path(source)
    if groups is None:
        path = source / 'structure/all_models.json'
        if not path.exists():
            raise ValueError('Previous sample lacks structure/all_models.json: ' + str(source))
        data = json.loads(path.read_text(encoding='utf-8'))
        if data.get('format') == 'orsearcher-model-groups-v1':
            # Shared evidence may represent many terminal variants. Shallow
            # expansion avoids duplicating every homology/domain payload;
            # subsequent evidence.apply replaces, rather than mutates, fields.
            groups = {cid: [dict(data['payloads'][p['_sequence_evidence_ref']],
                                 **{k: v for k, v in p.items() if k != '_sequence_evidence_ref'})
                            for p in placements] for cid, placements in data['groups'].items()}
        else:
            groups = data
    audit = {r['candidate_id']: r for r in table(source / 'annotation_qc/identity_audit.tsv')}
    selected = {r.get('selected_model_id') for r in audit.values()}
    reviewed = {r.get('selected_model_id') for r in audit.values() if r.get('structure_confidence') == 'reviewed_complete'}
    unique = {}
    for cid, entries in groups.items():
        for original in entries:
            m = {k: v for k, v in original.items() if k != 'discovery_alternatives'}
            m.setdefault('cid', cid)
            m['_baseline_selected'] = m['id'] in selected
            m['_reviewed_complete'] = m['id'] in reviewed
            m['_group_ids'] = {cid}
            m['_tools'] = {m.get('alignment_tool', m.get('source', 'unknown').split('_')[0])}
            if thresholds is not None:m['_thresholds']=dict(thresholds)
            m['_assessment'] = assessment(m, core_ids)
            key = placement(m)
            if key not in unique:
                unique[key] = m
            else:
                old = unique[key]
                merged_ids, merged_tools = old['_group_ids'] | m['_group_ids'], old['_tools'] | m['_tools']
                baseline, strict = old['_baseline_selected'] or m['_baseline_selected'], old['_reviewed_complete'] or m['_reviewed_complete']
                if rank(m) > rank(old):
                    unique[key] = m
                unique[key].update(_group_ids=merged_ids, _tools=merged_tools, _baseline_selected=baseline, _reviewed_complete=strict)
    pseudo = []
    eligible_pseudo_ids = {r['candidate_id'] for r in table(source / 'nonredundant/OR_catalogue.tsv')
                          if r.get('category') == 'pseudogene' and r.get('confidence') in ('genome_supported', 'assembly_supported_candidate')}
    for cid, decision in (decisions or {}).items():
        if decision.get('category') == 'pseudogene' and decision.get('confidence') == 'genome_supported':
            eligible_pseudo_ids.add(cid)
            selected.add(decision['model']['id'])
    for cid in eligible_pseudo_ids:
        entries = groups.get(cid, [])
        desired = (decisions or {}).get(cid, {}).get('model', {}).get('id') or audit.get(cid, {}).get('selected_model_id')
        found = next((m for m in entries if m['id'] == desired), None)
        if found is not None:
            m = dict(found)
            m.update(_group_ids={cid}, _tools={m.get('alignment_tool', 'unknown')})
            if thresholds is not None:m['_thresholds']=dict(thresholds)
            m['_assessment'] = assessment(m, core_ids)
            pseudo.append(m)
    return list(unique.values()), pseudo


def repair_plan(clusters, cross_status, refs, core_ids, args):
    planned, unavailable = [], []
    for cluster in clusters:
        m, a = cluster[0], cluster[0]['_assessment']
        group_ids = set().union(*(v['_group_ids'] for v in cluster))
        timedout = any(cross_status.get(cid, {}).get('status') in ('timeout', 'completed_with_warnings', 'error') for cid in group_ids)
        if a['usable_complete'] and not timedout:
            continue
        hits = {}
        for v in cluster:
            for h in v.get('homology_hits', []):
                if (h['ref'] in core_ids or m.get('_thresholds',{}).get('reconstruction_reference_policy')=='active-external') and h['ref'] in refs and h['bits'] >= tuning_policy(m)['repair_hit_bits'] and h['evalue'] <= tuning_policy(m)['repair_hit_evalue'] and h['identity'] >= tuning_policy(m)['repair_hit_identity']:
                    hits[h['ref']] = max(hits.get(h['ref'], 0), h['bits'])
        rids, seen = [], set()
        for rid in sorted(hits, key=lambda r: (-hits[r], r)):
            if refs[rid].sequence not in seen:
                rids.append(rid)
                seen.add(refs[rid].sequence)
            if len(rids) >= getattr(args, 'polish_references', 2):
                break
        cid = 'POLISH_' + digest(placement(m))[:20]
        if not rids:
            unavailable.append(dict(candidate_id=cid, source_model_id=m['id'], reason='no_complete_reconstruction_reference'))
            continue
        if bounds(m)[1] - bounds(m)[0] + 1 > getattr(args, 'max_window_bp', 500000):
            unavailable.append(dict(candidate_id=cid, source_model_id=m['id'], reason='model_span_over_window_limit'))
            continue
        priority = (not a['usable_complete'], a['reference_coverage'] >= .5, a['bits'], a['identity_pct'])
        planned.append(dict(candidate_id=cid, model=m, references=rids, priority=priority,
                            original_group_ids=sorted(group_ids), prior_timeout_or_error=timedout))
    planned.sort(key=lambda job: (job['priority'], job['candidate_id']), reverse=True)
    budget = getattr(args, 'polish_repair_budget', 256)
    return planned[:budget], planned[budget:], unavailable


def repair_case(job, root, args, refs, lengths):
    from .support import review_or_catalogue as review, crosscheck_or_catalogue as cross
    m, cid = job['model'], job['candidate_id']
    case = Path(root) / 'repair/cases' / cid
    case.mkdir(parents=True, exist_ok=True)
    chosen = {rid: refs[rid] for rid in job['references']}
    start, end = bounds(m)
    max_window = int(getattr(args, 'max_window_bp', 500000))
    span = end - start + 1
    if span > max_window:
        return [], dict(candidate_id=cid, status='unavailable', reason='model_span_over_window_limit')
    flank = min(int(getattr(args, 'polish_flank_bp', 80000)), max(0, (max_window - span) // 2))
    timeout = getattr(args, 'polish_timeout', 180)
    fingerprint = digest(dict(genome=metadata(args.genome), placement=placement(m),
                              references=[(rid, chosen[rid].sequence) for rid in sorted(chosen)],
                              flank_bp=flank, max_window_bp=max_window, max_intron=args.max_intron,
                              timeout=timeout, miniprot_timeout=min(60, timeout),
                              policy='practical_polish_repair_v1'))
    status_path, models_path = case / 'status.json', case / 'models.json'
    if status_path.exists() and models_path.exists():
        cached = json.loads(status_path.read_text())
        if cached.get('status') == 'completed' and cached.get('fingerprint') == fingerprint:
            values = json.loads(models_path.read_text())
            for new in values:
                new.update(_group_ids=set(job['original_group_ids']), _tools={new.get('alignment_tool', 'miniprot')})
            return values, dict(cached, cache_reused=True)
    row = dict(candidate_id=cid, sequence_id=m['contig'], strand=m['strand'], start=start, end=end,
               frameshifts=m.get('frameshifts', 0))
    parameters = SimpleNamespace(genome=args.genome, flank_bp=flank,
                                 samtools='samtools', miniprot='miniprot', threads_per_job=1,
                                 max_intron=args.max_intron, miniprot_timeout=min(60, timeout))
    values, issues = [], []
    try:
        values = review.model_candidates(row, case, [(s, e, 0) for s, e in m['exons']],
                                         m['protein'], chosen, lengths, parameters, original_model=m)
    except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
        issues.append(dict(tool='miniprot', status='timeout' if isinstance(error, subprocess.TimeoutExpired) else 'error', error=str(error)))
    region, region_metadata = case / 'region.fna', case / 'region.json'
    if region.exists() and region_metadata.exists():
        sequence = next(read_fasta(region)).sequence
        offset = json.loads(region_metadata.read_text())['offset']
        command = ['exonerate', '--model', 'protein2genome', '--querytype', 'protein', '--targettype', 'dna',
                   '--query', str(case / 'references.faa'), '--target', str(region), '--bestn', '3',
                   '--score', '100', '--proteinwordlen', '6', '--maxintron', str(args.max_intron),
                   '--refine', 'none', '--dpmemory', '128', '--showalignment', 'no',
                   '--showtargetgff', 'yes', '--showvulgar', 'yes', '--ryo', 'CDS_BEGIN\n%tcs\nCDS_END\n']
        save(case / 'exonerate_command.json', dict(command=command, timeout_seconds=timeout))
        output = case / 'exonerate.out'
        try:
            review.command(command, output, timeout)
            for alignment in review.parse_output(output):
                if alignment['strand'] != m['strand'] or not alignment['exons']:
                    continue
                ref = chosen[alignment['reference_id']]
                new = review.build_model(cid, 'polish_exonerate', sequence, offset, alignment['exons'], m['strand'],
                                         frameshifts=alignment['frameshifts'], ref=ref.id, qstart=alignment['qstart'] + 1,
                                         qend=alignment['qend'], ref_length=len(ref.sequence), identity=alignment['identity'])
                if new['cds'] != alignment['reported_cds'].upper() and not new['frameshifts']:
                    continue
                if review.overlap_fraction(m['exons'], new['exons']) < .45:
                    continue
                new.update(alignment_tool='exonerate', alignment=dict(qstart=alignment['qstart'] + 1,
                           qend=alignment['qend'], qlen=len(ref.sequence)))
                values.extend([new, *review.terminal_variants(new, sequence, offset, chosen)])
            save(case / 'exonerate_frame_events.json', cross.vulgar_frames(output, offset, sequence))
        except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
            issues.append(dict(tool='exonerate', status='timeout' if isinstance(error, subprocess.TimeoutExpired) else 'error', error=str(error)))
    kept = []
    for new in values:
        if new.get('source') == 'original' or overlap(m, new)[0] < .45:
            continue
        new.update(contig=m['contig'], cid=cid, annotation_gene_ids=[], fusion_pairs=[], bridge=False,
                   overlap=overlap(m, new)[0], _group_ids=set(job['original_group_ids']),
                   _tools={new.get('alignment_tool', 'miniprot')})
        kept.append(new)
    status = dict(candidate_id=cid, status='completed_with_warnings' if issues else 'completed',
                  references=list(chosen), model_count=len(kept), issues=issues, fingerprint=fingerprint,
                  region_flank_bp=flank, max_window_bp=max_window)
    save(case / 'models.json', [{k: v for k, v in model.items() if not k.startswith('_')} for model in kept])
    save(case / 'status.json', status)
    return kept, status


CATALOGUE_FIELDS = ['candidate_id', 'sequence_id', 'strand', 'start', 'end', 'category', 'confidence',
                    'independent_locus', 'overlapping_intact_model_ids', 'overlapping_pseudogene_model_ids', 'length_aa', 'length_cds',
                    'source_model_id', 'original_candidate_ids', 'coding_alternatives', 'best_reference',
                    'reference_coverage', 'query_coverage', 'identity_pct', 'alignment_tools',
                    'has_start', 'has_stop', 'frameshifts', 'internal_stops', 'identity_status',
                    'identity_reason', 'or_homology_bits', 'gr_homology_bits', 'or_profile_supported',
                    'gr_profile_supported', 'other_protein_bits', 'canonical_splicing',
                    'ambiguous_bases', 'strict_reviewed_complete', 'repair_derived', 'damage_evidence', 'experimentally_verified', 'completeness_blockers','completeness_policy','full_length_template_supported']


def export_set(folder, clusters, pseudo, args):
    from .engine import validate
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    prefix = re.sub(r'[^A-Za-z0-9_.-]', '_', getattr(args, 'model_prefix', None) or args.sample)
    entries = [(c[0], c, 'complete' if c[0]['_assessment']['usable_complete'] else 'partial') for c in clusters]
    seen = set()
    for m in pseudo:
        if placement(m) not in seen and m['_assessment']['identity_ok']:
            entries.append((m, [m], 'pseudogene'))
            seen.add(placement(m))
    entries.sort(key=lambda entry: (entry[0]['contig'], bounds(entry[0])[0], entry[0]['strand'], entry[2], entry[0]['id']))
    ids = {id(m): f'{prefix}_ORP{i:05d}' for i, (m, _, _) in enumerate(entries, 1)}
    pseudo_anchors, pseudo_relations = [], {}
    for m, _, category in entries:
        if category != 'pseudogene':
            continue
        intact_links = [ids[id(c[0])] for c in clusters if m['contig'] == c[0]['contig'] and m['strand'] == c[0]['strand'] and overlap(m, c[0])[2] >= 30]
        pseudo_links = [ids[id(a)] for a in pseudo_anchors if m['contig'] == a['contig'] and m['strand'] == a['strand'] and overlap(m, a)[2] >= 30]
        pseudo_relations[id(m)] = intact_links, pseudo_links
        if not intact_links and not pseudo_links:
            pseudo_anchors.append(m)
    rows, all_gff, sequences, cdss = [], ['##gff-version 3'], defaultdict(list), defaultdict(list)
    category_gff = {kind: ['##gff-version 3'] for kind in ('complete', 'partial', 'pseudogene')}
    for m, members, category in entries:
        cid, a = ids[id(m)], m['_assessment']
        start, end = bounds(m)
        overlaps, pseudo_overlaps = pseudo_relations.get(id(m), ([], []))
        independent = not bool(overlaps or pseudo_overlaps)
        confidence = ('usable_coding_supported' if a.get('completeness_policy')=='coding-and-identity' and not a.get('full_length_template_supported') else 'usable_external_homology_supported' if a.get('completeness_policy')=='external-homology' else 'usable_supported') if category == 'complete' else 'credible_partial' if category == 'partial' else 'assembly_supported_candidate'
        aliases = sorted(set().union(*(x.get('_group_ids', {x.get('cid', x['id'])}) for x in members)))
        signals=judge(m)['signals']
        row = dict(candidate_id=cid, sequence_id=m['contig'], strand=m['strand'], start=start, end=end,
                   category=category, confidence=confidence, independent_locus=independent,
                   overlapping_intact_model_ids=','.join(overlaps), overlapping_pseudogene_model_ids=','.join(pseudo_overlaps),
                   length_aa=m['length'], length_cds=len(m['cds']),
                   source_model_id=m['id'], original_candidate_ids=','.join(aliases), coding_alternatives=len(members) - 1,
                   best_reference=a['best_reference'], reference_coverage=round(a['reference_coverage'], 4),
                   query_coverage=round(a['query_coverage'], 4), identity_pct=a['identity_pct'],
                   alignment_tools=','.join(sorted(set().union(*(x.get('_tools', set()) for x in members)))),
                   has_start=bool(m.get('has_start')), has_stop=bool(m.get('has_stop')), frameshifts=m.get('frameshifts', 0),
                   internal_stops=len(m.get('stops', [])), identity_status=a['identity_status'], identity_reason=a['identity_reason'],
                   or_homology_bits=signals['or_homology_bits'],gr_homology_bits=signals['gr_homology_bits'],
                   or_profile_supported=signals['or_profile'],gr_profile_supported=signals['gr_profile'],
                   other_protein_bits=signals['background_bits'],canonical_splicing=bool(m.get('canonical')),
                   ambiguous_bases=m.get('ambiguous',0),strict_reviewed_complete=bool(m.get('_reviewed_complete')),
                   repair_derived=str(m.get('source','')).startswith('polish_') or str(m.get('cid','')).startswith('POLISH_'),
                   damage_evidence=('inherited_assembly_supported_not_revalidated_by_polish' if getattr(args,'command','')=='polish'
                                    else 'upstream_assembly_supported_structural_review') if category=='pseudogene' else 'not_applicable',
                   experimentally_verified=False, completeness_blockers=';'.join(a.get('completeness_blockers', [])),completeness_policy=a.get('completeness_policy','qualified-template'),full_length_template_supported=a.get('full_length_template_supported',False))
        rows.append(row)
        description = f'category={category} confidence={confidence} independent_locus={str(independent).lower()} damage_evidence={row["damage_evidence"]} source_model={m["id"]} predicted=true'
        sequences[category].append(FastaRecord(cid, description, m['protein']))
        cdss[category].append(FastaRecord(cid, description, m['cds']))
        attrs = f'ID={quote(cid, safe="_.-")};category={category};confidence={confidence}'
        gff = ['\t'.join(map(str, [m['contig'], 'ORsearcher', 'mRNA', start, end, '.', m['strand'], '.', attrs]))]
        n = 0
        for j, (s, e) in enumerate(sorted(m['exons'], reverse=m['strand'] == '-'), 1):
            gff.append('\t'.join(map(str, [m['contig'], 'ORsearcher', 'CDS', s, e, '.', m['strand'], (3 - n % 3) % 3,
                                               f'ID={quote(cid, safe="_.-")}.cds{j};Parent={quote(cid, safe="_.-")}'])))
            n += e - s + 1
        all_gff.extend(gff)
        category_gff[category].extend(gff)
    write_table(folder / 'OR_catalogue.tsv', rows, CATALOGUE_FIELDS)
    (folder / 'OR_models.gff3').write_text('\n'.join(all_gff) + '\n', encoding='utf-8')
    for category in ('complete', 'partial', 'pseudogene'):
        write_fasta(folder / f'OR_{category}.proteins.faa', sequences[category])
        write_fasta(folder / f'OR_{category}.cds.fna', cdss[category])
        (folder / f'OR_{category}.gff3').write_text('\n'.join(category_gff[category]) + '\n', encoding='utf-8')
    validation = validate(args.genome, folder)
    summary = dict(category_counts={kind: sum(r['category'] == kind for r in rows) for kind in ('complete', 'partial', 'pseudogene')},
                   independent_loci=sum(r['independent_locus'] for r in rows), non_independent_damage_alternatives=sum(not r['independent_locus'] for r in rows),
                   validation=validation, computational_prediction_not_experimental_truth=True)
    save(folder / 'summary.json', summary)
    return summary


def polish(args, sample_out, refs, evidence, groups=None, decisions=None):
    """Polish a completed annotation or a read-only previous sample directory."""
    sample_out = Path(sample_out)
    source = Path(getattr(args, 'previous_sample', None) or sample_out)
    folder, audit = sample_out / 'results', sample_out / 'polish_audit'
    audit.mkdir(parents=True, exist_ok=True)
    core_ids, core_policy = complete_reference_ids(args, refs)
    models, pseudo = load_pool(source, groups, decisions, core_ids, thresholds=tuning_policy(args=args))
    # Re-run cached evidence under the active reference fingerprint before
    # accepting any old label; replacing references cannot reuse stale labels.
    old_identity_policy = getattr(evidence.args, 'identity_policy', 'legacy')
    evidence.args.identity_policy = 'multichannel'
    try:
        all_models = models + pseudo
        prior_repeats = {id(model): (model.get('repeated_or_domains', False), model.get('repeated_reference_regions', False))
                         for model in all_models}
        for model in all_models:
            model.pop('background_hits', None)
        evidence.apply(all_models, 'polish_existing_' + args.sample)
        # Evidence.assign initializes repeat fields; a protein-score refresh
        # must not erase an earlier positive structural duplication warning.
        for model in all_models:
            old_domains, old_regions = prior_repeats[id(model)]
            model['repeated_or_domains'] = model.get('repeated_or_domains', False) or old_domains
            model['repeated_reference_regions'] = model.get('repeated_reference_regions', False) or old_regions
        background_evidence(all_models, args, audit / 'background_existing')
        for m in all_models:
            m['_assessment'] = assessment(m, core_ids)
        clusters, _, _ = consolidate(models)
        status_path = source / 'structure/crosscheck_summary.json'
        statuses = json.loads(status_path.read_text()).get('statuses', []) if status_path.exists() else []
        cross_status = {s['candidate_id']: s for s in statuses}
        selected, deferred, unavailable = repair_plan(clusters, cross_status, refs, core_ids, args)
        save(audit / 'repair_plan.json', dict(selected=[{k: v for k, v in job.items() if k != 'model'} for job in selected],
                                             deferred=[{k: v for k, v in job.items() if k != 'model'} for job in deferred], unavailable=unavailable,
                                             reference_policy=core_policy, complete_reference_ids=sorted(core_ids),
                                             not_a_whole_genome_search=True))
        LOG.info('Polish: coding models=%d initial loci=%d repair=%d deferred=%d', len(models), len(clusters), len(selected), len(deferred))
        extras, repair_statuses = [], []
        if selected:
            from .support import review_or_catalogue as review
            fai = Path(str(args.genome) + '.fai')
            if not fai.exists():
                review.command(['samtools', 'faidx', str(args.genome)])
            lengths = review.read_fai(fai)
            workers = max(1, min(int(getattr(args, 'workers', 1)), int(args.threads), len(selected)))
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(repair_case, job, audit, args, refs, lengths): job for job in selected}
                for index, future in enumerate(as_completed(futures), 1):
                    try:
                        values, status = future.result()
                    except Exception as error:
                        values, status = [], dict(candidate_id=futures[future]['candidate_id'], status='error', error=str(error))
                    extras.extend(values)
                    repair_statuses.append(status)
                    if index % 25 == 0 or index == len(selected):
                        LOG.info('Polish repair: %d/%d cases', index, len(selected))
            evidence.apply(extras, 'polish_repairs_' + args.sample)
            for model in extras:
                model.pop('background_hits', None)
            background_evidence(extras, args, audit / 'background_repairs')
            if extras:
                from .support import audit_or_loci as locus_audit
                repeat_args = SimpleNamespace(outdir=audit / 'repair', workers=int(args.threads),
                                              hmmsearch='hmmsearch', diamond='diamond', or_hmm=Path(args.or_hmm),
                                              or_db=args.or_db)
                try:
                    locus_audit.protein_repeats(extras, repeat_args)
                except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
                    # Missing a new fusion/repeat check is not permission to
                    # publish the new structure; old reviewed models survive.
                    for model in extras:
                        model['bridge'] = True
                    repair_statuses.append(dict(candidate_id='repair_structure_qc', status='error', error=str(error)))
        save(audit / 'repair_statuses.json', repair_statuses)
        combined = {placement(m): m for m in models}
        for m in extras:
            m['_assessment'] = assessment(m, core_ids)
            key = placement(m)
            if key in combined:
                combined[key]['_group_ids'].update(m['_group_ids'])
                combined[key]['_tools'].update(m['_tools'])
            else:
                combined[key] = m
        models = list(combined.values())
        clusters, assignments, conflicts = consolidate(models)
        accepted_pseudo = [m for m in pseudo if m['_assessment']['identity_ok']]
        write_table(audit / 'model_assignments.tsv', assignments)
        write_table(audit / 'competing_overlap_models.tsv', conflicts)
        chosen_ids = {c[0]['id'] for c in clusters} | {m['id'] for m in accepted_pseudo}
        alternate_rows = [dict(model_id=m['id'], contig=m['contig'], strand=m['strand'], start=bounds(m)[0], end=bounds(m)[1],
                               identity=m['_assessment']['identity_status'], reason=m['_assessment']['identity_reason'],
                               original_candidate_ids=','.join(sorted(m['_group_ids']))) for m in models if m['id'] not in chosen_ids]
        write_table(audit / 'alternative_models.tsv', alternate_rows,
                    ['model_id', 'contig', 'strand', 'start', 'end', 'identity', 'reason', 'original_candidate_ids'])
        # Retain physical coding paths plus union provenance, rather than
        # only a count. Repeated equivalent coordinate/CDS paths are deduped.
        from .model_store import save_groups
        saved_models = []
        for m in models + pseudo:
            saved = {k: v for k, v in m.items() if not k.startswith('_')}
            saved.update(polish_original_candidate_ids=sorted(m['_group_ids']), polish_alignment_tools=sorted(m['_tools']))
            saved_models.append(saved)
        save_groups(audit / 'models.json', {'polished_pool': saved_models})
        result = export_set(folder, clusters, accepted_pseudo, args)
        warnings = []
        if deferred:
            warnings.append('repair_budget_deferred')
        if unavailable:
            warnings.append('repair_window_or_reference_unavailable')
        if any(s['status'] != 'completed' for s in repair_statuses):
            warnings.append('local_repair_failures_retained')
        if conflicts:
            warnings.append('competing_overlap_paths_held_in_audit')
        if not getattr(args, 'background_reference', None):
            warnings.append('other_protein_background_not_checked')
        if not core_ids:
            warnings.append('no_complete_reconstruction_reference')
        if getattr(args,'command','')=='polish' and accepted_pseudo:
            warnings.append('pseudogene_damage_evidence_inherited_not_revalidated_by_polish')
        provenance=getattr(args,'previous_genome_provenance',None)
        if provenance is not None and not provenance['verified']:
            warnings.append('previous_genome_provenance_unverified')
        result.update(status='completed_with_warnings' if warnings else 'completed', warnings=warnings,
                      repair_selected=len(selected), repair_deferred=len(deferred), repair_unavailable=len(unavailable),
                      repair_status_counts=dict(Counter(s['status'] for s in repair_statuses)),
                      new_alignment_models=len(extras), alternative_models=len(alternate_rows),
                      source_sample=str(source), reference_completeness_policy=core_policy,
                      previous_genome_provenance=provenance,
                      policy='Practical multichannel OR identity and assembly-exact coding models; nontransitive locus consolidation; no reference-base grafting',
                      limitations=['Not an exhaustive OR census or experimentally verified gene set',
                                   'Overlapping competing paths may include real isoforms/overlapping genes and remain in the audit',
                                   'Assembly-supported pseudogenes may reflect assembly errors; overlapping intact alternatives are not independent loci',
                                   'Length and full-homology defaults can retain unusual ORs only in the candidate/audit sets'])
        save(folder / 'summary.json', result)
        save(audit / 'summary.json', result)
        LOG.info('Polish finished: complete=%d partial=%d pseudogene=%d independent=%d',
                 result['category_counts']['complete'], result['category_counts']['partial'],
                 result['category_counts']['pseudogene'], result['independent_loci'])
        return result
    finally:
        evidence.args.identity_policy = old_identity_policy
