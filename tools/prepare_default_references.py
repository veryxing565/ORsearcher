#!/usr/bin/env python3
"""Copy the public reference bundle or reproduce the legacy GR input locally."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "orsearcher" / "data" / "defaults"
SUPPLEMENT_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-73718-5/MediaObjects/41467_2026_73718_MOESM5_ESM.zip"
SUPPLEMENT_SHA256 = "909e06f7dafef8b67512796a5a7795b4445e70c7b203abc258ef6ae584970b7c"
RAW_GROR_SHA256 = "a702829d7cc44a05a0cf4035c2b0960167e6f887821a13f9dd2b796d454cd128"
LEGACY_GR_SHA256 = "80dc0466c4946d285220adedc10369caa20f220b989cd02cb97ab73efca3d02b"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def records(text: str):
    header = None
    seq = []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq).upper()
            header, seq = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError("Sequence appears before the first FASTA header")
            seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq).upper()


def extract_gr(archive: Path, selection: list[dict]) -> bytes:
    blob = archive.read_bytes()
    if sha256(blob) != SUPPLEMENT_SHA256:
        raise ValueError("Official supplement checksum differs from the pinned version")
    # Read one known member into memory; never extract arbitrary archive paths.
    with zipfile.ZipFile(archive) as handle:
        matches = [n for n in handle.namelist() if n.endswith("/sequences/unaligned_GROR.fasta")]
        if len(matches) != 1:
            raise ValueError("Cannot uniquely identify the original sequence FASTA")
        raw = handle.read(matches[0])
    if sha256(raw) != RAW_GROR_SHA256:
        raise ValueError("Original sequence FASTA checksum differs")
    source = {h.split()[0]: s.replace("?", "X") for h, s in records(raw.decode("utf-8-sig"))}
    output = []
    for item in selection:
        seq = source.get(item["source_id"])
        if seq is None or sha256(seq.encode("ascii")) != item["sequence_sha256"]:
            raise ValueError("Pinned GR sequence not reproduced: " + item["source_id"])
        output.append(">" + item["output_header"] + "\n")
        output.extend(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    return "".join(output).encode("ascii")


def build_selection(gr_file: Path, archive: Path, destination: Path):
    """Maintainer-only production of public ID/hash metadata, not sequences."""
    if destination.exists():
        raise FileExistsError("Refusing to replace existing selection metadata")
    with zipfile.ZipFile(archive) as handle:
        name = next(n for n in handle.namelist() if n.endswith("/sequences/unaligned_GROR.fasta"))
        raw = handle.read(name)
    if sha256(archive.read_bytes()) != SUPPLEMENT_SHA256 or sha256(raw) != RAW_GROR_SHA256:
        raise ValueError("Pinned official archive/FASTA checksum mismatch")
    source = {h.split()[0]: s.replace("?", "X") for h, s in records(raw.decode("utf-8-sig"))}
    selection = []
    for header, sequence in records(gr_file.read_text(encoding="utf-8")):
        fields = dict(part.split("=", 1) for part in header.split()[1:] if "=" in part)
        original = fields["original_id"]
        if source.get(original) != sequence:
            raise ValueError("Selected reference differs from original normalized sequence: " + original)
        selection.append(dict(source_id=original, output_header=header, sequence_sha256=sha256(sequence.encode("ascii"))))
    if len(selection) != 412 or extract_gr(archive, selection) != gr_file.read_bytes():
        raise ValueError("The exact 412-sequence GR panel was not reproduced")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(selection, indent=2) + "\n", encoding="utf-8")


def install_bytes(destination: Path, data: bytes, expected: str):
    if sha256(data) != expected:
        raise ValueError("Pinned reference checksum mismatch: " + destination.name)
    if destination.exists():
        if sha256(destination.read_bytes()) != expected:
            raise FileExistsError("Refusing to replace a different existing file: " + str(destination))
        return
    temporary = destination.with_name(destination.name + ".partial")
    temporary.write_bytes(data)
    temporary.replace(destination)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, help="Destination for a copy of the public bundle or the legacy GR input")
    parser.add_argument("--gr-file", type=Path, help="Import the exact previously prepared 412-sequence GR.faa")
    parser.add_argument("--gr-archive", type=Path, help="Reproduce the local subset from the pinned original official ZIP")
    parser.add_argument("--download-gr", action="store_true", help="Download the original official supplement, then reproduce the local subset")
    parser.add_argument("--build-selection", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--accept-gr-license", action="store_true", help="Acknowledge CC BY-NC-ND 4.0 and take responsibility for permitted local use; does not grant redistribution rights")
    options = parser.parse_args(argv)
    if options.build_selection:
        if not (options.gr_file and options.gr_archive and options.accept_gr_license):
            parser.error("Maintainer metadata build requires --gr-file, --gr-archive and --accept-gr-license")
        build_selection(options.gr_file, options.gr_archive, options.build_selection)
        print("Pinned selection metadata written; exact local GR reconstruction verified")
        return 0
    if not options.outdir:
        parser.error("--outdir is required")
    if sum(bool(x) for x in (options.gr_file, options.gr_archive, options.download_gr)) > 1:
        parser.error("Choose only one of --gr-file, --gr-archive and --download-gr")
    manifest = json.loads((DATA / "manifest.json").read_text(encoding="utf-8"))
    options.outdir.mkdir(parents=True, exist_ok=True)
    gr_source = options.gr_file or options.gr_archive or options.download_gr
    if gr_source and not options.accept_gr_license:
        parser.error("Read docs/DATA_PROVENANCE.md, then explicitly pass --accept-gr-license for the restricted GR resource")
    try:
        if not gr_source:
            for name, item in manifest['files'].items():
                install_bytes(options.outdir / name, (DATA / name).read_bytes(), item['sha256'])
            for name in ('manifest.json','DATA_LICENSES.md'):
                blob = (DATA / name).read_bytes()
                install_bytes(options.outdir / name, blob, sha256(blob))
            if 'or_sources' in manifest:
                source = manifest['or_sources']
                name = source['file']
                install_bytes(options.outdir / name, (DATA / name).read_bytes(), source['sha256'])
            print('Verified public bundle:', options.outdir.resolve())
            return 0
        if options.gr_file:
            gr_blob = options.gr_file.read_bytes()
        else:
            archive = options.gr_archive
            if options.download_gr:
                archive = options.outdir / 'official_mitchell_2026.zip'
                if not archive.exists():
                    request = urllib.request.Request(SUPPLEMENT_URL, headers={'User-Agent':'ORsearcher/2.3-reference-preparation'})
                    with urllib.request.urlopen(request, timeout=90) as response:
                        install_bytes(archive, response.read(), SUPPLEMENT_SHA256)
            selection = json.loads((DATA / 'GR_selection.json').read_text(encoding='utf-8'))
            gr_blob = extract_gr(archive, selection)
        target = options.outdir / 'legacy_Mitchell_GR.faa'
        install_bytes(target, gr_blob, LEGACY_GR_SHA256)
        print('Verified legacy GR input:', target.resolve())
        print('Use --gr-reference with the legacy input and the other original run settings.')
        return 0
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        print("Reference preparation failed: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
