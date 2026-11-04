# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Validate Exonerate coding coordinates and export same-locus evidence."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from .audit_rescue_endpoints import table, transcript, endpoint_options, write_table
from .rescue_flagged_or import load_candidate_exons, overlap_fraction
from orsearcher.fasta import FastaRecord, read_fasta, translate, write_fasta


def parse_output(path):
    records, current = [], None
    reading_cds = False
    in_dump = False
    for line in path.read_text().splitlines():
        if line.startswith("vulgar: "):
            fields = line.split()
            current = dict(reference_id=fields[1], qstart=int(fields[2]), qend=int(fields[3]),
                           strand=fields[8], score=int(fields[9]), exons=[], reported_cds="",
                           frameshifts=sum(fields[i] == "F" for i in range(10, len(fields), 3)), identity=0.)
            records.append(current)
            reading_cds = False
        elif line == "# --- START OF GFF DUMP ---":
            in_dump = True
        elif line == "# --- END OF GFF DUMP ---":
            in_dump = False
        elif current is not None and in_dump and not line.startswith("#"):
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "cds":
                current["exons"].append((int(fields[3]), int(fields[4])))
            elif len(fields) == 9 and fields[2] == "gene":
                match = re.search(r"(?:^|;)\s*identity\s+([0-9.]+)", fields[8])
                if match:
                    current["identity"] = float(match.group(1)) / 100
        elif current is not None and line == "CDS_BEGIN":
            reading_cds = True
        elif line == "CDS_END":
            reading_cds = False
        elif current is not None and reading_cds:
            current["reported_cds"] += line.strip()
    return records


