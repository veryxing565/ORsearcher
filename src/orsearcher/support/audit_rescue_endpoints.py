# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Audit genomic disruptions and propose bounded terminal-exon extensions.

Proposals are not experimentally confirmed gene models. No genomic bases are
changed, internal stops are retained, and originals are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from orsearcher.fasta import FastaRecord, read_fasta, reverse_complement, translate, write_fasta
from orsearcher.models import attach_reference_lengths, parse_miniprot_gff
from orsearcher.sharding import Segment

STOPS = {"TAA", "TAG", "TGA"}


def table(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def transcript(sequence, exons, strand, leading_phase=0):
    positions = []
    for start, end in sorted(exons, reverse=strand == "-"):
        positions.extend(range(start - 1, end) if strand == "+" else range(end - 1, start - 2, -1))
    positions = positions[leading_phase:]
    bases = "".join(sequence[position] for position in positions)
    if strand == "-":
        # Positions already run 5' -> 3'; complement without reversing again.
        bases = reverse_complement(bases)[::-1]
    return bases, positions


def oriented_bases(sequence, positions, strand):
    if not positions or min(positions) < 0 or max(positions) >= len(sequence):
        return ""
    bases = "".join(sequence[position] for position in positions)
    return reverse_complement(bases)[::-1] if strand == "-" else bases


def endpoint_options(sequence, cds, positions, strand, max_extension=450, max_trim=90):
    """Return start adjustments in nt and the first in-frame terminal stop.

    Negative start adjustment extends upstream. Positive trims inside the first
    coding exon. Upstream searches stop at an in-frame stop or ambiguous base.
    """
    direction = 1 if strand == "+" else -1
    starts = []
    if cds.startswith("ATG"):
        starts.append(0)
    else:
        for shift in range(-3, -max_extension - 1, -3):
            codon_positions = [positions[0] + direction * (shift + i) for i in range(3)]
            codon = oriented_bases(sequence, codon_positions, strand)
            if len(codon) != 3 or any(base not in "ACGT" for base in codon) or codon in STOPS:
                break
            if codon == "ATG":
                starts.append(shift)
        for shift in range(3, min(max_trim, len(cds) - 3) + 1, 3):
            # Trimming cannot jump an intron or an upstream termination codon.
            if positions[shift] != positions[0] + direction * shift:
                break
            codon = cds[shift:shift + 3]
            if any(base not in "ACGT" for base in codon) or codon in STOPS:
                break
            if codon == "ATG":
                starts.append(shift)
    if len(cds) % 3 == 0 and cds[-3:] in STOPS:
        return starts, 0
    for added in range(1, max_extension + 1):
        if (len(cds) + added) % 3:
            continue
        added_positions = [positions[-1] + direction * i for i in range(1, added + 1)]
        bases = oriented_bases(sequence, added_positions, strand)
        if len(bases) != added or any(base not in "ACGT" for base in bases):
            return starts, None
        if (cds + bases)[-3:] in STOPS:
            return starts, added
    return starts, None


