"""Reference FASTA and HMMER3 inputs for one analysis."""
from dataclasses import dataclass
from pathlib import Path
import gzip
import math

from .fasta import FastaRecord, detect_sequence_type, read_fasta, write_fasta
from .runtime import save


PROTEIN = frozenset('ACDEFGHIKLMNPQRSTVWYBXZJUO')


@dataclass(frozen=True)
class Profile:
    name: str
    accession: str
    text: str
    ga: tuple[float, float] | None


def read_profiles(path):
    """Read one or more amino-acid HMMER3 profiles by content."""
    path = Path(path)
    try:
        raw = path.read_bytes()
        if raw.startswith(b'\x1f\x8b'):
            raw = gzip.decompress(raw)
        lines = raw.decode('ascii').splitlines()
    except (OSError, ValueError, EOFError) as error:
        raise ValueError(f'Cannot read HMMER3 profiles from {path}: {error}') from error
    profiles = []
    block = []
    for line in lines:
        if not block and not line.strip():
            continue
        block.append(line)
        if line.strip() != '//':
            continue
        fields = {}
        for row in block[1:-1]:
            parts = row.split(maxsplit=1)
            if not parts:
                continue
            key = parts[0]
            if key == 'HMM':
                break
            fields[key] = parts[1].strip() if len(parts) > 1 else ''
        try:
            name = fields['NAME']
            if (not block[0].startswith('HMMER3/') or not name or len(name.split()) != 1 or
                    fields['ALPH'].lower() != 'amino' or int(fields['LENG']) <= 0 or
                    not any(row.split()[:1] == ['HMM'] for row in block[1:])):
                raise ValueError('invalid amino-acid profile header')
            ga = tuple(map(float, fields['GA'].rstrip(';').split())) if 'GA' in fields else None
            if ga is not None and (len(ga) != 2 or not all(map(math.isfinite, ga))):
                raise ValueError('invalid GA cutoffs')
        except (KeyError, ValueError) as error:
            raise ValueError(f'{path}: expected a valid amino-acid HMMER3 profile ({error})') from error
        profiles.append(Profile(name, fields.get('ACC', ''), '\n'.join(block)+'\n', ga))
        block = []
    if block or not profiles:
        raise ValueError(f'{path}: expected complete HMMER3 profiles ending with //')
    return profiles


def merge_profiles(paths):
    profiles = {}
    for path in paths:
        for model in read_profiles(path):
            existing = profiles.get(model.name)
            if existing is not None and existing.text != model.text:
                raise ValueError(f'{path}: conflicting HMM NAME {model.name!r}. Use unique profile names within a family.')
            profiles[model.name] = model
    return list(profiles.values())


def user_proteins(path):
    for rec in read_fasta(path):
        seq = rec.sequence.removesuffix('*')
        if not seq or set(seq) - PROTEIN:
            raise ValueError(f'{path}: {rec.id!r} must contain amino acids without gaps or internal stops')
        if len(seq) < 50:
            raise ValueError(f'{path}: {rec.id!r} is shorter than 50 amino acids')
        if detect_sequence_type(seq) == 'dna':
            raise ValueError(f'{path}: {rec.id!r} looks like DNA. Supply translated protein FASTA.')
        yield FastaRecord(rec.id, rec.description, seq)


def merge_proteins(family, base, replacement, additions):
    if replacement is None:
        records = list(read_fasta(base))
        sources = additions
    else:
        records = []
        sources = [replacement, *additions]
    seen = {r.sequence.rstrip('*') for r in records}
    used = {r.id for r in records}
    mapping = []
    for source_index, path in enumerate(sources, 1):
        for record_index, rec in enumerate(user_proteins(path), 1):
            if rec.sequence in seen:
                mapping.append(dict(file=str(path), original_header=rec.description, status='duplicate_sequence'))
                continue
            rid = f'user_{family}_{source_index:02d}_{record_index:06d}'
            while rid in used:
                rid += '_'
            records.append(FastaRecord(rid, rec.description, rec.sequence))
            mapping.append(dict(file=str(path), original_header=rec.description, reference_id=rid, status='used'))
            seen.add(rec.sequence)
            used.add(rid)
    return records, mapping


def materialize_references(args, root):
    """Keep normalized inputs in the output tree, leaving original inputs untouched."""
    folder = root/'shared/reference_inputs'
    folder.mkdir(parents=True, exist_ok=True)
    report = {}
    args.user_or_ids = set()
    for family in ('or', 'gr'):
        field = f'{family}_reference'
        path = getattr(args, field)
        replacement = getattr(args, f'replace_{family}', None)
        additions = getattr(args, f'add_{family}', [])
        if replacement or additions:
            records, mapping = merge_proteins(family, path, replacement, additions)
            active = folder/f'{family}.faa'
            if not active.exists():
                write_fasta(active, records)
            setattr(args, field, active)
            report[family+'_proteins'] = mapping
            if family == 'or':
                args.user_or_ids = {m['reference_id'] for m in mapping if m['status'] == 'used'}
        hmm_paths = getattr(args, 'profile_inputs', {}).get(family)
        if hmm_paths:
            models = merge_profiles(hmm_paths)
            active = folder/f'{family}.hmm'
            if not active.exists():
                active.write_text(''.join(p.text for p in models), encoding='ascii')
            setattr(args, f'{family}_hmm', active if family=='or' else [active])
            report[family+'_profiles'] = [dict(name=p.name, accession=p.accession) for p in models]
    save(folder/'selection.json', dict(sources=args.reference_selection, normalized=report))
