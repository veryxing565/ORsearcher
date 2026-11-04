"""Annotation comparison, separate from prediction and biological validation.

Coordinates are 1-based inclusive. Names select the comparator, not predictions.
An unmatched prediction is not automatically a false positive.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path

from .models import parse_attributes
from .runtime import save
from .support.review_or_catalogue import write_table

BENCHMARK_SCHEMA = 2
GENE_TYPES = {'gene', 'pseudogene'}
TRANSCRIPT_TYPES = {'mRNA', 'transcript', 'pseudogenic_transcript'}


def receptor_family(name):
    """Conservative symbols/products; arbitrary identifiers need a family map."""
    name = str(name).strip()
    stem = r'(Orco|Or\d+[a-z]*(?:[_.-]\d+)?|Gr\d+[a-z]*(?:[_.-]\d+)?)'
    match = re.fullmatch(stem, name, re.I)
    if not match:
        match = re.fullmatch(r'[A-Za-z][A-Za-z0-9]{1,12}[_-]' + stem, name, re.I)
    if not match:
        # Require a capitalized family boundary for unseparated species prefixes.
        match = re.fullmatch(r'[A-Za-z][A-Za-z0-9]{1,9}(Orco|ORCO|(?:OR|Or|GR|Gr)\d+[A-Za-z]*(?:[_.-]\d+)?)', name)
    if match:
        return 'GR' if match.group(1).lower().startswith('gr') else 'OR'
    product = name.replace('_', ' ')
    if re.search(r'\b(?:binding|associated|accessory|trafficking|regulator)\b', product, re.I):
        return ''
    families = set()
    if re.search(r'\b(?:odorant|olfactory)[ -]+receptor\b', product, re.I):
        families.add('OR')
    if re.search(r'\bgustatory[ -]+receptor\b', product, re.I):
        families.add('GR')
    return next(iter(families)) if len(families) == 1 else ''


def gff_lines(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            yield line


def _features(path):
    for number, line in enumerate(gff_lines(path), 1):
        if not line.strip() or line.startswith('#'):
            continue
        fields = line.rstrip('\r\n').split('\t')
        if len(fields) != 9:
            raise ValueError(f'{path}:{number}: expected nine GFF3 columns')
        try:
            start, end = int(fields[3]), int(fields[4])
        except ValueError as error:
            raise ValueError(f'{path}:{number}: invalid coordinates') from error
        if start < 1 or end < start:
            raise ValueError(f'{path}:{number}: invalid coordinate interval')
        attrs = parse_attributes(fields[8])
        if fields[2] in GENE_TYPES | TRANSCRIPT_TYPES and not attrs.get('ID'):
            raise ValueError(f'{path}:{number}: gene/transcript requires GFF3 ID')
        if fields[2] == 'CDS' and not attrs.get('Parent'):
            raise ValueError(f'{path}:{number}: CDS requires GFF3 Parent (GTF is not supported)')
        yield line if line.endswith('\n') else line + '\n', fields, attrs


def _aliases(attrs):
    values = set()
    for key in ('ID', 'gene', 'Name', 'gene_name', 'locus_tag', 'protein_id', 'transcript_id', 'Alias', 'Dbxref'):
        values.update(x.strip() for x in attrs.get(key, '').split(',') if x.strip())
    return values


def _family_map(path):
    if not path:
        return {}
    if isinstance(path, dict):
        return path
    result = {}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'feature_id', 'family'} <= set(reader.fieldnames or []):
            raise ValueError('Family map requires TSV columns: feature_id, family, optional name')
        for row in reader:
            key, family = row['feature_id'].strip(), row['family'].strip().upper()
            if not key or family not in ('OR', 'GR', 'OTHER'):
                raise ValueError('Family map requires nonempty feature_id and OR/GR/OTHER family')
            item = dict(family=family, name=row.get('name', '').strip())
            if key in result and result[key] != item:
                raise ValueError('Conflicting family-map rows: ' + key)
            result[key] = item
    return result


def _normalize_exons(exons):
    """Union for coverage only; retain CDS boundaries/phases for structure."""
    union = []
    for start, end in sorted(set(exons)):
        if union and start <= union[-1][1] + 1:
            union[-1] = (union[-1][0], max(end, union[-1][1]))
        else:
            union.append((start, end))
    return union


def _coverage(reference, prediction):
    reference, prediction = _normalize_exons(reference), _normalize_exons(prediction)
    shared = sum(max(0, min(b, d) - max(a, c) + 1) for a, b in reference for c, d in prediction)
    return (shared / sum(b - a + 1 for a, b in reference),
            shared / sum(d - c + 1 for c, d in prediction))


def annotation_models(path, family_map=None, or_name_patterns=(), gr_name_patterns=()):
    """Select named or explicitly mapped coding loci, including CDS products.

    Explicit maps override automatic names for their locus. Unknown, ambiguous
    or conflicting map entries fail instead of changing the denominator silently.
    Stream unrelated features instead of keeping the whole annotation in memory.
    """
    mappings = _family_map(family_map)
    patterns = {'OR': [re.compile(p) for p in or_name_patterns],
                'GR': [re.compile(p) for p in gr_name_patterns]}
    nodes, parents, declared_genes = {}, {}, set()
    node_parts = defaultdict(list)
    for _, fields, attrs in _features(path):
        if fields[2] not in GENE_TYPES | TRANSCRIPT_TYPES:
            continue
        key = attrs['ID']
        if key in nodes:
            old_fields, old_attrs = nodes[key]
            multipart = ('part' in attrs or 'part' in old_attrs or
                         attrs.get('exception') == 'trans-splicing' or old_attrs.get('exception') == 'trans-splicing')
            if not multipart or old_fields[2] != fields[2] or old_attrs.get('Parent', '') != attrs.get('Parent', ''):
                raise ValueError('Duplicate gene/transcript ID without compatible multipart annotation: ' + key)
            for label in ('gene', 'Name'):
                if old_attrs.get(label) != attrs.get(label):
                    raise ValueError('Inconsistent multipart annotation labels: ' + key)
        else:
            nodes[key] = (fields, attrs)
            parents[key] = [p for p in attrs.get('Parent', '').split(',') if p]
        node_parts[key].append((fields, attrs))
        if fields[2] in GENE_TYPES:
            declared_genes.add(key)
    ancestry = {}

    def gene_of(key, trail=()):
        if key in ancestry:
            return ancestry[key]
        if key in trail:
            raise ValueError('Cyclic annotation Parent chain: ' + key)
        if key in declared_genes:
            ancestry[key] = key
            return key
        if key not in nodes:
            raise ValueError('Annotation CDS/transcript Parent not found: ' + key)
        gene_ids = {gene_of(parent, trail + (key,)) for parent in parents[key]}
        if len(gene_ids) != 1:
            raise ValueError('Transcript must resolve to one gene: ' + key)
        ancestry[key] = next(iter(gene_ids))
        return ancestry[key]

    transcripts = {key: gene_of(key) for key in nodes if key not in declared_genes}
    hints, explicit, used_map = defaultdict(set), defaultdict(list), defaultdict(set)
    coding, model_info, partials = defaultdict(set), {}, defaultdict(bool)
    for _, fields, attrs in _features(path):
        kind = fields[2]
        if kind not in GENE_TYPES | TRANSCRIPT_TYPES | {'CDS'}:
            continue
        tids = ([p for p in attrs.get('Parent', '').split(',') if p]
                if kind == 'CDS' else [attrs['ID']])
        gids = {gene_of(tid) for tid in tids}
        aliases = _aliases(attrs)
        auto = {receptor_family(value) for value in aliases}
        auto.update(receptor_family(attrs.get(key, '')) for key in ('product', 'description'))
        for family, regexes in patterns.items():
            if any(pattern.search(value) for pattern in regexes for value in aliases):
                auto.add(family)
        auto.discard('')
        mapped = [(key, mappings[key]) for key in sorted(aliases & mappings.keys())]
        for gid in gids:
            hints[gid].update(auto)
            for key, entry in mapped:
                if isinstance(entry, str):
                    entry = dict(family=entry, name='')
                explicit[gid].append(entry)
                used_map[key].add(gid)
            if kind in GENE_TYPES and (attrs.get('partial', '').lower() == 'true' or
                                      attrs.get('pseudo', '').lower() == 'true' or kind == 'pseudogene' or
                                      'start_range' in attrs or 'end_range' in attrs):
                partials[gid] = True
        if kind != 'CDS':
            continue
        for tid in tids:
            gid = gene_of(tid)
            start, end = int(fields[3]), int(fields[4])
            phase = fields[7]
            if phase not in ('0', '1', '2', '.'):
                raise ValueError('Invalid CDS phase: ' + tid)
            coding[tid].add((start, end, phase))
            partial = attrs.get('partial', '').lower() == 'true'
            partial |= attrs.get('pseudo', '').lower() == 'true'
            partial |= 'start_range' in attrs or 'end_range' in attrs
            if tid in nodes:
                partial |= nodes[tid][1].get('partial', '').lower() == 'true'
                partial |= nodes[tid][1].get('pseudo', '').lower() == 'true'
                partial |= 'start_range' in nodes[tid][1] or 'end_range' in nodes[tid][1]
            info = model_info.setdefault(tid, dict(partial=False, locations=set()))
            info['partial'] |= partial
            info['locations'].add((fields[0], fields[6], start, end))
    unused = set(mappings) - used_map.keys()
    ambiguous = [key for key, gids in used_map.items() if len(gids) != 1]
    if unused:
        raise ValueError('Family-map IDs absent from annotation: ' + ', '.join(sorted(unused)[:10]))
    if ambiguous:
        raise ValueError('Family-map IDs match multiple genes: ' + ', '.join(sorted(ambiguous)[:10]))
    genes = {}
    for gid in sorted(declared_genes):
        families = {entry['family'] for entry in explicit[gid]} if explicit[gid] else hints[gid]
        if len(families) > 1:
            raise ValueError('Conflicting OR/GR family labels for gene: ' + gid)
        family = next(iter(families), '')
        if family not in ('OR', 'GR'):
            continue
        fields, attrs = nodes[gid]
        positions = {(part[0][0], part[0][6]) for part in node_parts[gid]}
        if len(positions) != 1 or fields[6] not in ('+', '-'):
            raise ValueError('Selected receptor has unsupported multi-contig/strand gene parts: ' + gid)
        names = [entry.get('name', '') for entry in explicit[gid] if entry.get('name')]
        genes[gid] = dict(gene_id=gid, name=names[0] if names else attrs.get('gene', attrs.get('Name', gid)),
                          family=family, contig=fields[0], strand=fields[6],
                          annotation_partial=str(partials[gid]).lower(),
                          family_source='explicit_map' if explicit[gid] else 'annotation_name_or_product')
    models = []
    for tid, parts in sorted(coding.items()):
        gid = gene_of(tid)
        if gid not in genes:
            continue
        ordered = sorted(parts)
        gene_parts = node_parts[gid]
        lower = min(int(part[0][3]) for part in gene_parts)
        upper = max(int(part[0][4]) for part in gene_parts)
        for seqid, strand, start, end in model_info[tid]['locations']:
            if seqid != genes[gid]['contig'] or strand != genes[gid]['strand']:
                raise ValueError('Annotation CDS contig/strand inconsistent with gene: ' + tid)
            if start < lower or end > upper:
                raise ValueError('Annotation CDS outside gene bounds: ' + tid)
        if len({(start, end) for start, end, _ in ordered}) != len(ordered):
            raise ValueError('Conflicting phases at identical annotation CDS coordinates: ' + tid)
        models.append(dict(genes[gid], transcript_id=tid,
                           exons=[(start, end) for start, end, _ in ordered],
                           phases=[phase for _, _, phase in ordered],
                           annotation_partial=str(partials[gid] or model_info[tid]['partial'] or
                               (ordered[-1 if genes[gid]['strand'] == '-' else 0][2] in ('1', '2'))).lower()))
    selected = {key for key in nodes if gene_of(key) in genes}
    features = []
    for line, fields, attrs in _features(path):
        if attrs.get('ID') in selected or any(parent in selected for parent in attrs.get('Parent', '').split(',')):
            features.append((line, fields, attrs))
    return genes, transcripts, models, features


def _annotation(args):
    return annotation_models(args.annotation, getattr(args, 'family_map', None),
                             getattr(args, 'or_name_pattern', ()) or (),
                             getattr(args, 'gr_name_pattern', ()) or ())


def prepare(args):
    genes, _, models, features = _annotation(args)
    if not any(m['family'] == 'OR' for m in models):
        raise ValueError('No coding OR comparator recognized; supply --family-map or naming patterns')
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'receptor_reference.gff3').write_text('##gff-version 3\n' + ''.join(row[0] for row in features), encoding='utf-8')
    save(out / 'annotation_inventory.json', dict(genes=dict(Counter(g['family'] for g in genes.values())),
         coding_transcripts=dict(Counter(m['family'] for m in models)), annotation=str(args.annotation),
         limitation='Annotation comparison only; not independent functional validation'))
    return out / 'receptor_reference.gff3'


def _prediction_folder(args):
    root = Path(args.prediction)
    sample = getattr(args, 'sample', None)
    if sample:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', sample) or sample in ('.', '..'):
            raise ValueError('Invalid benchmark --sample')
        root = root / 'samples' / sample if (root / 'samples').is_dir() else root / sample
    elif (root / 'samples').is_dir():
        samples = sorted(p for p in (root / 'samples').iterdir() if p.is_dir())
        if len(samples) != 1:
            raise ValueError('Batch prediction requires --sample (multiple or zero samples)')
        root = samples[0]
    scope = getattr(args, 'prediction_scope', 'auto')
    if scope not in ('auto', 'final', 'legacy', 'candidates'):
        raise ValueError('Unknown prediction scope: ' + str(scope))
    options = []
    if scope in ('auto', 'final'):
        options += [(root / 'results', 'final', 'OR_catalogue.tsv', 'OR_models.gff3')]
    if scope in ('auto', 'final', 'legacy'):
        options += [(root, 'final' if root.name == 'results' or scope == 'final' else 'provided_models', 'OR_catalogue.tsv', 'OR_models.gff3')]
    if scope in ('auto', 'legacy'):
        options += [(root / 'nonredundant', 'legacy', 'OR_catalogue.tsv', 'OR_models.gff3')]
    if scope in ('auto', 'candidates'):
        options += [(root, 'candidates', 'candidate_catalogue.tsv', 'candidate_models.gff3')]
    for folder, kind, catalogue, gff in options:
        present = [(folder / name).is_file() for name in (catalogue, gff)]
        if any(present):
            if not all(present):
                raise ValueError('Incomplete prediction catalogue/GFF pair: ' + str(folder))
            return folder, kind, folder / catalogue, folder / gff
    raise FileNotFoundError('No compatible prediction pair found under ' + str(root))


def _seqid_map(path):
    if not path:
        return {}
    result = {}
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'prediction_seqid', 'annotation_seqid'} <= set(reader.fieldnames or []):
            raise ValueError('Seqid map requires prediction_seqid and annotation_seqid columns')
        for row in reader:
            key, value = row['prediction_seqid'].strip(), row['annotation_seqid'].strip()
            if not key or not value or (key in result and result[key] != value):
                raise ValueError('Empty/conflicting seqid mapping')
            result[key] = value
    return result


def _prediction_models(catalogue_path, gff_path, seqids):
    with catalogue_path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'candidate_id', 'sequence_id', 'strand', 'category', 'confidence'} <= set(reader.fieldnames or []):
            raise ValueError('Prediction catalogue lacks required fields')
        rows = list(reader)
    catalogue = {}
    for row in rows:
        cid = row['candidate_id']
        if not cid or cid in catalogue:
            raise ValueError('Empty/duplicate candidate ID: ' + cid)
        if row['strand'] not in ('+', '-'):
            raise ValueError('Invalid prediction strand: ' + cid)
        value = row.get('independent_locus', 'true').lower()
        if value not in ('true', 'false', '1', '0'):
            raise ValueError('Invalid independent_locus flag: ' + cid)
        row['_independent'] = value in ('true', '1')
        catalogue[cid] = row
    parts = defaultdict(set)
    for _, fields, attrs in _features(gff_path):
        if fields[2] != 'CDS':
            continue
        for cid in attrs['Parent'].split(','):
            if cid not in catalogue:
                raise ValueError('GFF CDS absent from catalogue: ' + cid)
            row = catalogue[cid]
            if fields[0] != row['sequence_id'] or fields[6] != row['strand']:
                raise ValueError('Prediction catalogue/GFF contig or strand mismatch: ' + cid)
            if fields[7] not in ('0', '1', '2', '.'):
                raise ValueError('Invalid prediction CDS phase: ' + cid)
            parts[cid].add((int(fields[3]), int(fields[4]), fields[7]))
    for cid, row in catalogue.items():
        if not parts[cid]:
            raise ValueError('Prediction missing CDS: ' + cid)
        ordered = sorted(parts[cid])
        if len({(a, b) for a, b, _ in ordered}) != len(ordered):
            raise ValueError('Conflicting prediction CDS phases: ' + cid)
        row['_exons'] = [(a, b) for a, b, _ in ordered]
        row['_phases'] = [phase for _, _, phase in ordered]
        for key, boundary in (('start', ordered[0][0]), ('end', max(b for _, b, _ in ordered))):
            if row.get(key) not in (None, ''):
                try:
                    consistent = int(row[key]) == boundary
                except ValueError:
                    consistent = False
                if not consistent:
                    raise ValueError('Prediction catalogue/GFF bounds mismatch: ' + cid)
        row['_seqid'] = seqids.get(row['sequence_id'], row['sequence_id'])
        row['_overlapping_cds'] = any(right[0] <= left[1] for left, right in zip(ordered, ordered[1:]))
    unknown = set(seqids) - {r['sequence_id'] for r in rows}
    if unknown:
        raise ValueError('Seqid-map prediction IDs absent from catalogue: ' + ', '.join(sorted(unknown)))
    return rows


def _matching(edges, minimum, support, exact=False, complete=False):
    """Maximum-cardinality matching: duplicates/fusions cannot inflate recall."""
    graph = defaultdict(dict)
    for hit in edges:
        if hit['reference_family'] != 'OR' or not hit['independent_locus']:
            continue
        if hit['reference_CDS_covered'] < minimum or hit['predicted_CDS_supported'] < support:
            continue
        if exact and not hit['splice_and_endpoints_match']:
            continue
        if complete and (hit['category'] != 'complete' or hit['annotation_partial'] == 'true'):
            continue
        cid, gid = hit['candidate_id'], hit['gene_id']
        score = (hit['reference_CDS_covered'] * hit['predicted_CDS_supported'], hit['transcript_id'])
        previous = graph[cid].get(gid)
        if previous is None or score > previous[0]:
            graph[cid][gid] = (score, hit)
    assigned = {}

    def augment(cid, visited):
        for gid in sorted(graph[cid], key=lambda key: (-graph[cid][key][0][0], key)):
            if gid in visited:
                continue
            visited.add(gid)
            if gid not in assigned or augment(assigned[gid], visited):
                assigned[gid] = cid
                return True
        return False

    for cid in sorted(graph):
        augment(cid, set())
    return {gid: graph[cid][gid][1] for gid, cid in sorted(assigned.items())}


def _sha256(path):
    if not path:
        return None
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def evaluate(args):
    minimum = getattr(args, 'min_reference_coverage', .5)
    support = getattr(args, 'min_predicted_support', .3)
    if not 0 < minimum <= 1 or not 0 < support <= 1:
        raise ValueError('Benchmark coverage/support thresholds must be in (0,1]')
    genes, _, gold, _ = _annotation(args)
    coding_ids = {model['gene_id'] for model in gold}
    ors = {key: gene for key, gene in genes.items() if gene['family'] == 'OR' and key in coding_ids}
    if not ors:
        raise ValueError('No coding OR comparator recognized; supply --family-map or naming patterns')
    folder, scope, catalogue_path, gff_path = _prediction_folder(args)
    rows = _prediction_models(catalogue_path, gff_path, _seqid_map(getattr(args, 'seqid_map', None)))
    candidate_pool = scope == 'candidates'
    byregion = defaultdict(list)
    for model in gold:
        byregion[model['contig'], model['strand']].append(model)
    # Non-receptor contigs are still valid comparison space. A prediction there
    # is unmatched, not evidence of an assembly/identifier mismatch.
    reference_seqids = {fields[0] for _, fields, _ in _features(args.annotation)}
    prediction_seqids = {row['_seqid'] for row in rows}
    if rows and not reference_seqids & prediction_seqids:
        raise ValueError('No shared prediction/annotation seqids; check assembly and explicit --seqid-map')
    edges, hits = [], []
    qualifying, reference_predictions = defaultdict(set), defaultdict(set)
    for row in rows:
        pe, matches = row['_exons'], []
        for model in byregion[row['_seqid'], row['strand']]:
            a, b = _coverage(model['exons'], pe)
            if not a or not b:
                continue
            ge = model['exons']
            gj = {(x[1] + 1, y[0] - 1) for x, y in zip(ge, ge[1:])}
            pj = {(x[1] + 1, y[0] - 1) for x, y in zip(pe, pe[1:])}
            near = gj == pj and abs(pe[0][0] - ge[0][0]) <= 3 and abs(pe[-1][1] - ge[-1][1]) <= 3
            phase_match = row['_phases'] == model['phases'] and '.' not in row['_phases']
            exact = pe == ge and phase_match and not row['_overlapping_cds']
            hit = dict(candidate_id=row['candidate_id'], category=row['category'], confidence=row['confidence'],
                       independent_locus=row['_independent'], gene_id=model['gene_id'], gene_name=model['name'],
                       reference_family=model['family'], transcript_id=model['transcript_id'],
                       annotation_partial=model['annotation_partial'], reference_CDS_covered=a, predicted_CDS_supported=b,
                       CDS_phase_match=int(phase_match), splice_endpoints_within_3bp=int(near),
                       splice_and_endpoints_match=int(exact))
            edges.append(hit)
            matches.append(hit)
            if model['family'] == 'OR' and a >= minimum and b >= support and row['_independent']:
                qualifying[row['candidate_id']].add(model['gene_id'])
                reference_predictions[model['gene_id']].add(row['candidate_id'])
        if matches:
            best = sorted(matches, key=lambda h: (-h['reference_CDS_covered'] * h['predicted_CDS_supported'], h['gene_id'], h['transcript_id']))[0]
            hits.append(dict(best))
        else:
            hits.append(dict(candidate_id=row['candidate_id'], category=row['category'], confidence=row['confidence'],
                             independent_locus=row['_independent'], reference_family='unmatched'))
    covered = _matching(edges, minimum, support)
    tiers = {percentage: _matching(edges, percentage / 100, support) for percentage in (50, 75, 90)}
    exact = _matching(edges, minimum, support, exact=True, complete=True) if not candidate_pool else {}
    candidate_exact = _matching(edges, minimum, support, exact=True)
    complete_coordinate_exact = _matching([h for h in edges if h['category'] == 'complete'], minimum, support, exact=True)
    gr_overlap = {hit['candidate_id'] for hit in edges if hit['reference_family'] == 'GR'
                  and hit['reference_CDS_covered'] >= minimum and hit['predicted_CDS_supported'] >= support
                  and hit['independent_locus']}
    for hit in hits:
        cid = hit['candidate_id']
        hit['qualifying_OR_loci'] = len(qualifying[cid])
        hit['assigned_OR_gene_id'] = next((gid for gid, h in covered.items() if h['candidate_id'] == cid), '')
        hit['also_overlaps_GR'] = cid in gr_overlap
    out = Path(args.outdir)
    if out.resolve() == folder.resolve() or folder.resolve() in out.resolve().parents or out.resolve() in folder.resolve().parents:
        raise ValueError('Benchmark output must be separate from the prediction tree')
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise FileExistsError('Benchmark output must be empty: ' + str(out))
    fields = ['candidate_id', 'category', 'confidence', 'independent_locus', 'gene_id', 'gene_name', 'reference_family',
              'transcript_id', 'annotation_partial', 'reference_CDS_covered', 'predicted_CDS_supported',
              'CDS_phase_match', 'splice_endpoints_within_3bp', 'splice_and_endpoints_match']
    write_table(out / 'candidate_annotation_comparison.tsv', hits,
                fields + ['qualifying_OR_loci', 'assigned_OR_gene_id', 'also_overlaps_GR'])
    write_table(out / 'all_annotation_overlaps.tsv', edges, fields)
    missing = [dict(ors[key], reason='no_one_to_one_candidate_above_declared_coverage_support') for key in sorted(ors) if key not in covered]
    write_table(out / 'missed_OR_annotation.tsv', missing,
                ['gene_id', 'name', 'family', 'contig', 'strand', 'annotation_partial', 'family_source', 'reason'])
    unscorable = [dict(gene, reason='no_annotated_CDS_model') for gid, gene in sorted(genes.items()) if gid not in coding_ids]
    write_table(out / 'unscorable_reference_loci.tsv', unscorable,
                ['gene_id', 'name', 'family', 'contig', 'strand', 'annotation_partial', 'family_source', 'reason'])
    locus_rows = []
    for gid in sorted(ors):
        assigned = covered.get(gid, {})
        locus_rows.append(dict(gene_id=gid, gene_name=ors[gid]['name'],
                               assigned_candidate=assigned.get('candidate_id', ''),
                               reference_CDS_covered=assigned.get('reference_CDS_covered', ''),
                               recovered_50pct=gid in tiers[50], recovered_75pct=gid in tiers[75], recovered_90pct=gid in tiers[90],
                               complete_exact_structure=gid in exact,
                               qualifying_prediction_records=len(reference_predictions[gid])))
    write_table(out / 'OR_locus_recovery.tsv', locus_rows,
                ['gene_id', 'gene_name', 'assigned_candidate', 'reference_CDS_covered', 'recovered_50pct', 'recovered_75pct',
                 'recovered_90pct', 'complete_exact_structure', 'qualifying_prediction_records'])
    summary = dict(benchmark_schema_version=BENCHMARK_SCHEMA, prediction_scope='unverified_discovery_candidates' if candidate_pool else scope,
                   resolved_prediction_directory=str(folder.resolve()), reference_OR_gene_loci=len(ors),
                   reference_OR_annotated_gene_loci=sum(g['family'] == 'OR' for g in genes.values()),
                   reference_OR_loci_without_CDS=sum(g['family'] == 'OR' for g in unscorable),
                   reference_OR_coding_transcripts=sum(m['family'] == 'OR' for m in gold),
                   predicted_representative_records=len(rows), predicted_independent_records=sum(r['_independent'] for r in rows),
                   OR_annotation_loci_recovered=len(covered), locus_recall_against_annotation=len(covered) / len(ors),
                   OR_loci_with_complete_splice_endpoint_match=len(exact), predictions_overlapping_GR=len(gr_overlap),
                   OR_loci_with_coordinate_phase_match=len(candidate_exact),
                   OR_loci_with_complete_prediction_coordinate_phase_match=len(complete_coordinate_exact),
                   reference_GR_coding_loci=len({m['gene_id'] for m in gold if m['family'] == 'GR'}),
                   annotation_unmatched_predictions=sum(h['reference_family'] == 'unmatched' for h in hits),
                   merged_OR_prediction_records=sum(len(gids) > 1 for gids in qualifying.values()),
                   duplicated_OR_reference_loci=sum(len(cids) > 1 for cids in reference_predictions.values()),
                   prediction_records_with_overlapping_CDS=sum(r['_overlapping_cds'] for r in rows),
                   prediction_seqids_absent_from_comparator=sorted(prediction_seqids - reference_seqids),
                   confidence_counts=dict(Counter(r['confidence'] for r in rows)),
                   parameters=dict(min_reference_coverage=minimum, min_predicted_support=support,
                                   or_name_patterns=getattr(args, 'or_name_pattern', ()) or (),
                                   gr_name_patterns=getattr(args, 'gr_name_pattern', ()) or (),
                                   matching='deterministic_maximum_cardinality_one_to_one', exact_endpoint_tolerance_bp=0,
                                   exact_requires_CDS_phase=True),
                   input_sha256={key: _sha256(path) for key, path in dict(annotation=args.annotation,
                       prediction_catalogue=catalogue_path, prediction_gff=gff_path,
                       family_map=getattr(args, 'family_map', None), seqid_map=getattr(args, 'seqid_map', None)).items()},
                   limitations=['Annotation is a comparator, not absolute biological truth',
                       'Unmatched predictions are not automatically false positives; no precision/FDR is inferred',
                       'Annotation selection by names/products does not establish comparator completeness; inventory must be audited',
                       'Coordinate/phase agreement is not sequence, functional or biological full-length validation',
                       'A predicted fusion can claim at most one reference locus; all overlaps are retained',
                       'Coverage thresholds are independently matched; locus assignments may differ between tiers',
                       'Reference partial/pseudogene loci do not count as intact complete structure recovery',
                       'No target annotation supplied by this scorer to prediction; reference leakage must be controlled separately'])
    for percentage, matching in tiers.items():
        summary[f'OR_loci_recovered_at_{percentage}pct'] = len(matching)
    if not summary['reference_GR_coding_loci']:
        summary['limitations'].append('No coding GR comparator recognized; zero GR overlaps is not evidence of GR specificity')
    if candidate_pool:
        summary.update(records_are_not_gene_counts=True,
                       OR_loci_with_candidate_splice_endpoint_match=len(candidate_exact))
        summary['limitations'].append('Recall-only candidates do not confirm OR identity, full length or pseudogene status')
    save(out / 'benchmark_summary.json', summary)
    return summary
