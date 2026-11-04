#!/usr/bin/env python3
"""Verify the packaged reference files and their release manifest."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import tarfile
import zipfile


DATA = '/orsearcher/data/defaults/'
REQUIRED = {'OR.faa', 'GR.faa', 'complete_OR.faa', 'OR.hmm', 'GR_1.hmm', 'GR_2.hmm', 'background.faa'}
LEGACY_GR_SHA256 = '80dc0466c4946d285220adedc10369caa20f220b989cd02cb97ab73efca3d02b'


def verify(path):
    archive = zipfile.ZipFile(path) if zipfile.is_zipfile(path) else tarfile.open(path)
    with archive:
        if isinstance(archive,zipfile.ZipFile):
            names = archive.namelist()
            read = archive.read
        else:
            names = archive.getnames()
            def read(name):
                with archive.extractfile(name) as handle:
                    return handle.read()
        manifests = [name for name in names if ('/'+name).endswith(DATA+'manifest.json')]
        assert len(manifests)==1, (path, manifests)
        manifest = json.loads(read(manifests[0]))
        directory = manifests[0].rsplit('/',1)[0]
        assert set(manifest['files'])==REQUIRED and manifest['schema']=='orsearcher-reference-panel-v2'
        assert manifest['files']['GR.faa']['license']=='CC-BY-4.0'
        assert manifest['files']['GR.faa']['sha256']!=LEGACY_GR_SHA256
        for name, info in manifest['files'].items():
            data = read(f'{directory}/{name}')
            assert len(data)==info['bytes'] and sha256(data).hexdigest()==info['sha256'], (path, name)
        if 'or_sources' in manifest:
            info = manifest['or_sources']
            table = read(f"{directory}/{info['file']}")
            assert sha256(table).hexdigest() == info['sha256']
            assert len(table.splitlines())-1 == manifest['files']['OR.faa']['records']
        for name in names:
            if name.lower().endswith(('.fa','.fna','.faa','.fasta')):
                assert sha256(read(name)).hexdigest()!=LEGACY_GR_SHA256, (path, name)
        assert f'{directory}/DATA_LICENSES.md' in names
    print(f'{path.name}: verified {len(REQUIRED)} reference files and public GR provenance')


if __name__=='__main__':
    for item in sys.argv[1:]:
        verify(Path(item))
