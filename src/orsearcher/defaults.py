"""Select the baseline reference inputs and run-local additions."""
from __future__ import annotations

import json
import os
from pathlib import Path

FILES = {
    'or_reference': 'OR.faa', 'gr_reference': 'GR.faa',
    'or_hmm': 'OR.hmm', 'gr_hmm': ['GR_1.hmm', 'GR_2.hmm'],
    'background_reference': 'background.faa', 'complete_reference': 'complete_OR.faa',
}


def reference_directory(value=None):
    return Path(value or os.environ.get('ORSEARCHER_REFERENCE_DIR') or
                Path(__file__).parent / 'data' / 'defaults').expanduser().resolve()


def inventory(value=None):
    root = reference_directory(value)
    result = {'directory': str(root), 'files': {}}
    for field, names in FILES.items():
        for name in names if isinstance(names, list) else [names]:
            path = root / name
            result['files'][name] = {'exists': path.is_file(),
                                      'bytes': path.stat().st_size if path.is_file() else 0}
    manifest = root / 'manifest.json'
    if manifest.is_file():
        result['provenance'] = json.loads(manifest.read_text(encoding='utf-8'))
    return result


def resolve_references(args):
    from .input_files import merge_profiles, user_proteins

    root = reference_directory(getattr(args, 'reference_dir', None))
    custom_or = any(getattr(args, name, None) is not None for name in ('or_reference','replace_or'))
    args.use_added_templates = getattr(args, 'complete_reference', None) is None
    if getattr(args, 'replace_or', None) is not None and not args.use_added_templates:
        raise ValueError('--replace-or infers template eligibility automatically. For an explicit subset, use --or-reference with --complete-reference.')
    changes = {}
    for flag in ('replace_hmm', 'add_hmm'):
        by_role = {'or': [], 'gr': []}
        for role, filename in getattr(args, flag, []) or []:
            if role.lower() not in by_role:
                raise ValueError(f'--{flag.replace("_", "-")} expects OR or GR before the filename')
            by_role[role.lower()].append(Path(filename))
        changes[flag] = by_role
    selected = {}
    args.profile_inputs = {}
    args.extra_reference_inputs = []
    for field, names in FILES.items():
        value = getattr(args, field, None)
        # A custom OR set must not inherit the unrelated bundled full-length set.
        if field == 'complete_reference' and value is None and (custom_or or getattr(args,'infer_templates',False)):
            selected[field] = 'audited_reconstruction_eligible_subset'
            continue
        if field == 'background_reference' and getattr(args, 'no_background', False):
            setattr(args, field, None)
            selected[field] = 'disabled'
            continue
        family = field.split('_',1)[0]
        replacement = getattr(args, 'replace_'+family, None) if field in ('or_reference','gr_reference') else None
        additions = getattr(args, 'add_'+family, []) if field in ('or_reference','gr_reference') else []
        if replacement is not None:
            if value is not None:
                raise ValueError(f'Use either --replace-{family} or --{family}-reference')
            value = replacement
        if field in ('or_hmm','gr_hmm') and changes['replace_hmm'][family]:
            if value is not None:
                raise ValueError(f'Use either --replace-hmm {family.upper()} or --{field.replace("_", "-")}')
            paths = changes['replace_hmm'][family]
            value = paths[0] if family=='or' else paths
        if value is None:
            value = [root / n for n in names] if isinstance(names, list) else root / names
            selected[field] = 'default'
        else:
            selected[field] = 'explicit_replacement'
        paths = value if isinstance(value, list) else [value]
        if field in ('or_hmm','gr_hmm'):
            repl = changes['replace_hmm'][family]
            extra = changes['add_hmm'][family]
            if repl or extra:
                args.profile_inputs[family] = (repl or paths)+extra
                merge_profiles(args.profile_inputs[family])
                args.extra_reference_inputs.extend(args.profile_inputs[family])
                if extra:
                    selected[field] += '_with_additions'
        for path in [*paths, *additions]:
            if not Path(path).is_file() or not Path(path).stat().st_size:
                raise ValueError(f'Missing/empty {field.replace("_", "-")}: {path}. '
                                 'Configure --reference-dir (or ORSEARCHER_REFERENCE_DIR), '
                                 'or supply the matching reference option; see docs/DATA_PROVENANCE.md.')
        for path in ([replacement] if replacement else [])+additions:
            for _ in user_proteins(path):
                pass
        if additions:
            args.extra_reference_inputs.extend(additions)
            selected[field] += '_with_additions'
        setattr(args, field, value)
    args.reference_selection = selected
    return selected
