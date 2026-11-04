"""Create a clean source ZIP with the verified public reference bundle.

Run from any directory. Does not push to GitHub or change a Git repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import zipfile


TOP_FILES = {'.gitignore', 'README.md', 'LICENSE', 'CITATION.cff', 'CONTRIBUTING.md',
             'CHANGELOG.md', 'environment.yml', 'pyproject.toml', 'MANIFEST.in'}
TOP_DIRS = {'src', 'tests', 'docs', 'tools', 'scripts', '.github', 'examples'}
DATA_DIR = 'src/orsearcher/data/defaults'


def release_files(root, bundled):
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        rel = path.relative_to(root)
        if rel.parts[0] not in TOP_DIRS and rel.as_posix() not in TOP_FILES:
            continue
        if any(p in {'__pycache__', '.git', '.pytest_cache', '.reference_build'} or p.endswith('.egg-info')
               or p.startswith('.reference_test_') for p in rel.parts):
            continue
        if rel.parent.as_posix() == DATA_DIR and path.name not in bundled:
            continue
        if path.suffix in {'.pyc', '.pyo', '.log', '.partial'}:
            continue
        yield path


def build(root, destination):
    root, destination = Path(root).resolve(), Path(destination).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f'Refusing to overwrite {destination}')
    manifest = json.loads((root / DATA_DIR / 'manifest.json').read_text(encoding='utf-8'))
    if manifest['files']['GR.faa']['license'] != 'CC-BY-4.0':
        raise ValueError('Expected the public UniProtKB GR panel')
    for name, info in manifest['files'].items():
        if hashlib.sha256((root / DATA_DIR / name).read_bytes()).hexdigest() != info['sha256']:
            raise ValueError('Reference checksum mismatch: ' + name)
    bundled = set(manifest['files']) | {'manifest.json','DATA_LICENSES.md','GR_selection.json'}
    if 'or_sources' in manifest:
        source = manifest['or_sources']
        if hashlib.sha256((root / DATA_DIR / source['file']).read_bytes()).hexdigest() != source['sha256']:
            raise ValueError('OR source table checksum mismatch')
        bundled.add(source['file'])
    files = list(release_files(root,bundled))
    with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as archive:
        for path in files:
            archive.write(path, 'orsearcher-2.3.0/' + path.relative_to(root).as_posix())
    report = {'archive': destination.name, 'files': len(files),
              'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
              'bundled_panel': manifest['panel_id']}
    destination.with_suffix('.manifest.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    build(Path(__file__).resolve().parents[1], args.out)
