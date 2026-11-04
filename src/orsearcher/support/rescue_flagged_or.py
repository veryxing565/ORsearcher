# Reused tested v1 coordinate/structure primitives; v2 orchestration is separate.
#!/usr/bin/env python3
"""Targeted, non-destructive genome re-alignment of flagged short OR models.

This is a post-processing diagnostic, not a replacement for the original
ORsearcher result. It uses full-length reference proteins and only searches
small windows around the candidate loci.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

from orsearcher.fasta import FastaRecord, read_fasta, write_fasta
from orsearcher.models import attach_reference_lengths, extract_sequences, parse_attributes, parse_miniprot_gff
from orsearcher.sharding import Segment


LOG = logging.getLogger("rescue_flagged_or")
FIELDS = (
    "candidate_id", "sequence_id", "original_length_aa", "reference_count",
    "alignment_count", "same_locus_count", "best_reference", "best_length_aa",
    "best_reference_coverage", "best_cds_overlap_fraction", "best_has_start",
    "best_has_stop", "best_frameshifts", "best_inframe_stops",
    "best_internal_stops", "best_or_hmm_ga", "best_gr_hmm_ga",
    "length_gain_aa", "decision", "notes",
)


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str], stdout_path: Path | None = None,
        stderr_path: Path | None = None) -> str:
    LOG.info("Command: %s", " ".join(command))
    kwargs: dict[str, object] = {"text": True, "check": False}
    if stdout_path:
        with stdout_path.open("w", encoding="utf-8") as stdout:
            result = subprocess.run(command, stdout=stdout, stderr=subprocess.PIPE, **kwargs)
    else:
        result = subprocess.run(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, **kwargs)
    if stderr_path:
        stderr_path.write_text(result.stderr or "", encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Command failed with exit {result.returncode}: {' '.join(command)}\n{result.stderr}")
    return result.stdout or ""


def read_fai(path: Path) -> dict[str, int]:
    contigs = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            contigs[fields[0]] = int(fields[1])
    return contigs


def load_candidate_exons(gff: Path, wanted: set[str]) -> dict[str, list[tuple[int, int]]]:
    exons = {name: [] for name in wanted}
    with gff.open(encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            parent = parse_attributes(fields[8]).get("Parent")
            if parent in exons:
                exons[parent].append((int(fields[3]), int(fields[4])))
    missing = [name for name, intervals in exons.items() if not intervals]
    if missing:
        raise ValueError(f"Candidates absent from input GFF CDS: {', '.join(missing)}")
    return exons


def overlap_fraction(original: list[tuple[int, int]], alternate: list[tuple[int, int]]) -> float:
    bases = sum(end - start + 1 for start, end in original)
    shared = sum(max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
                 for left_start, left_end in original
                 for right_start, right_end in alternate)
    return min(1.0, shared / bases) if bases else 0.0


def positive_hmm_ids(hmmsearch: str, model: Path, proteins: Path, prefix: Path) -> set[str]:
    table = prefix.with_suffix(".tbl")
    run([hmmsearch, "--cut_ga", "--noali", "--tblout", str(table),
         str(model), str(proteins)], stderr_path=prefix.with_suffix(".stderr.log"))
    return {line.split()[0] for line in table.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")}


def choose_references(hits: list[dict[str, str]], full_references: dict[str, FastaRecord],
                      query_id: str, limit: int) -> list[FastaRecord]:
    names = []
    for row in hits:
        if row["qseqid"] == query_id and row["sseqid"] in full_references:
            if row["sseqid"] not in names:
                names.append(row["sseqid"])
        if len(names) >= limit:
            break
    return [full_references[name] for name in names]


def process_one(row: dict[str, str], sequence: str, references: list[FastaRecord],
                exons: list[tuple[int, int]], contig_length: int,
                args: argparse.Namespace) -> dict[str, object]:
    candidate_id = row["candidate_id"]
    case_dir = args.outdir / candidate_id
    case_dir.mkdir(parents=True, exist_ok=True)
    start = max(1, int(row["start"]) - args.flank_bp)
    end = min(contig_length, int(row["end"]) + args.flank_bp)
    region = run([args.samtools, "faidx", str(args.genome),
                  f"{row['sequence_id']}:{start}-{end}"],
                 stderr_path=case_dir / "samtools.stderr.log")
    region_sequence = "".join(region.splitlines()[1:]).upper()
    if len(region_sequence) != end - start + 1:
        raise ValueError(f"Wrong extracted region length for {candidate_id}")
    write_fasta(case_dir / "region.fna", [FastaRecord("region", "", region_sequence)])
    write_fasta(case_dir / "full_length_references.faa", references)

    alignment = case_dir / "miniprot.gff3"
    run([args.miniprot, "-t", str(args.threads), "--gff", "-G", str(args.max_intron),
         "--outc", "0.15", "--outs", "0.20", "-N", "30", "--outn", "30",
         "-n", "3", "-p", "0.2", str(case_dir / "region.fna"),
         str(case_dir / "full_length_references.faa")],
        stdout_path=alignment, stderr_path=case_dir / "miniprot.stderr.log")

    segment = Segment("region", row["sequence_id"], "", start - 1, end,
                      contig_length, "region.fna", False)
    models = parse_miniprot_gff(alignment, {"region": segment})
    attach_reference_lengths(models, {record.id: len(record.sequence) for record in references})
    valid = []
    for model in models:
        try:
            extract_sequences(model, region_sequence)
        except Exception as error:
            LOG.warning("Skipping untranslatable model for %s: %s", candidate_id, error)
            continue
        if not model.protein_sequence or model.strand != row["strand"]:
            continue
        alternate_exons = [(start - 1 + feature.start, start - 1 + feature.end)
                           for feature in model.features if feature.feature_type == "CDS"]
        shared = overlap_fraction(exons, alternate_exons)
        if shared < args.min_cds_overlap:
            continue
        valid.append((model, shared))

    # HMMER cannot score literal stop codons; preserve the disruption count in
    # the TSV while representing those positions as unknown amino acids here.
    all_proteins = [FastaRecord(f"alt{index:04d}", "",
                                model.protein_sequence.replace("*", "X"))
                    for index, (model, _) in enumerate(valid, 1)]
    if all_proteins:
        proteins_path = case_dir / "same_locus_alternatives.faa"
        write_fasta(proteins_path, all_proteins)
        or_positive = positive_hmm_ids(args.hmmsearch, args.or_hmm, proteins_path,
                                        case_dir / "or_hmm")
        gr_positive = set()
        for index, model in enumerate(args.gr_hmm, 1):
            gr_positive.update(positive_hmm_ids(args.hmmsearch, model, proteins_path,
                                                 case_dir / f"gr_hmm_{index}"))
    else:
        or_positive, gr_positive = set(), set()

    details = []
    for index, (model, shared) in enumerate(valid, 1):
        alt_id = f"alt{index:04d}"
        details.append({
            "alt_id": alt_id, "reference_id": model.reference_id,
            "protein_length": len(model.protein_sequence),
            "reference_coverage": round(model.coverage, 4),
            "cds_overlap_fraction": round(shared, 4),
            "has_start": int(model.has_start), "has_stop": int(model.has_stop),
            "frameshifts": model.frameshifts,
            "inframe_stops": model.inframe_stops,
            "internal_stops": model.internal_stops,
            "or_hmm_ga": int(alt_id in or_positive),
            "gr_hmm_ga": int(alt_id in gr_positive),
            "genomic_start": model.original_start, "genomic_end": model.original_end,
            "exon_count": model.exon_count,
        })
    with (case_dir / "same_locus_alternatives.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(details[0]) if details else
                                ["alt_id", "reference_id", "protein_length"], delimiter="\t")
        writer.writeheader()
        writer.writerows(details)

    base = {name: "" for name in FIELDS}
    base.update(candidate_id=candidate_id, sequence_id=row["sequence_id"],
                original_length_aa=len(sequence), reference_count=len(references),
                alignment_count=len(models), same_locus_count=len(valid))
    if not details:
        base["decision"] = "no_same_locus_alignment"
        base["notes"] = "Keep original OR-family candidate; full length remains unverified"
        return base

    details.sort(key=lambda item: (
        item["or_hmm_ga"],
        not (item["frameshifts"] or item["inframe_stops"] or item["internal_stops"]),
        item["protein_length"], item["reference_coverage"], item["cds_overlap_fraction"]
    ), reverse=True)
    best = details[0]
    original_length = len(sequence)
    intact = not (best["frameshifts"] or best["inframe_stops"] or best["internal_stops"])
    full_length = (
        best["protein_length"] >= max(args.min_full_length_aa, int(0.8 *
                                                  next(len(ref.sequence) for ref in references
                                                       if ref.id == best["reference_id"])))
        and best["reference_coverage"] >= 0.75
        and best["has_start"] and best["has_stop"] and intact and best["or_hmm_ga"]
    )
    if full_length and best["protein_length"] >= original_length + args.min_gain_aa:
        decision = "provisional_full_length_rescue"
    elif best["protein_length"] >= original_length + args.min_gain_aa:
        decision = "extended_model_needs_review"
    else:
        decision = "no_meaningful_extension"
    base.update(
        best_reference=best["reference_id"], best_length_aa=best["protein_length"],
        best_reference_coverage=best["reference_coverage"],
        best_cds_overlap_fraction=best["cds_overlap_fraction"],
        best_has_start=best["has_start"], best_has_stop=best["has_stop"],
        best_frameshifts=best["frameshifts"],
        best_inframe_stops=best["inframe_stops"],
        best_internal_stops=best["internal_stops"],
        best_or_hmm_ga=best["or_hmm_ga"], best_gr_hmm_ga=best["gr_hmm_ga"],
        length_gain_aa=best["protein_length"] - original_length,
        decision=decision,
        notes="Diagnostic model only; original annotation was not overwritten",
    )
    return base


