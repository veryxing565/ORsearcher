from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable
from urllib.parse import unquote

from .fasta import reverse_complement, translate
from .sharding import Segment


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            attributes[unquote(key)] = unquote(value)
    return attributes


@dataclass
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str]


@dataclass
class Candidate:
    raw_id: str
    virtual_id: str
    original_id: str
    original_start: int
    original_end: int
    strand: str
    reference_id: str
    reference_start: int
    reference_end: int
    reference_length: int
    identity: float
    positive: float
    frameshifts: int
    inframe_stops: int
    rank: int
    shard: str
    edge_distance: int
    features: list[Feature] = field(default_factory=list)
    reference_hits: set[str] = field(default_factory=set)
    cds_sequence: str = ""
    protein_sequence: str = ""
    locus_sequence: str = ""
    has_start: bool = False
    has_stop: bool = False
    internal_stops: int = 0
    tm_helices: int = 0
    status: str = ""
    confidence: str = ""

    @property
    def coverage(self) -> float:
        return (self.reference_end - self.reference_start + 1) / max(1, self.reference_length)

    @property
    def span(self) -> int:
        return self.original_end - self.original_start + 1

    @property
    def exon_count(self) -> int:
        return sum(feature.feature_type == "CDS" for feature in self.features)

    def quality_key(self) -> tuple[float, ...]:
        intact = 1.0 if self.frameshifts == 0 and self.inframe_stops == 0 else 0.0
        return (
            intact,
            min(self.coverage, 1.0),
            self.identity,
            min(self.edge_distance, 100_000) / 100_000,
            -float(self.rank),
        )


TARGET_RE = re.compile(r"^(\S+)\s+(\d+)\s+(\d+)(?:\s+([+-]))?$")


def parse_miniprot_gff(path: Path, mapping: dict[str, Segment]) -> list[Candidate]:
    candidates: dict[str, Candidate] = {}
    with path.open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_no}: expected 9 GFF3 columns")
            seqid, source, feature_type, start, end, score, strand, phase, attr_text = fields
            if seqid not in mapping:
                raise ValueError(f"{path}:{line_no}: unknown shard sequence {seqid!r}")
            attrs = parse_attributes(attr_text)
            feature = Feature(
                seqid, source, feature_type, int(start), int(end), score, strand, phase, attrs
            )
            segment = mapping[seqid]
            if feature_type == "mRNA":
                raw_id = attrs.get("ID")
                target = attrs.get("Target", "")
                match = TARGET_RE.match(target)
                if not raw_id or not match:
                    raise ValueError(f"{path}:{line_no}: mRNA lacks usable ID/Target attributes")
                reference_id, qstart, qend, _ = match.groups()
                qstart_i, qend_i = int(qstart), int(qend)
                # Target length is absent in GFF. Query end is a conservative
                # denominator until the normalized reference length is attached.
                original_start = segment.original_start + feature.start
                original_end = segment.original_start + feature.end
                candidates[raw_id] = Candidate(
                    raw_id=raw_id,
                    virtual_id=seqid,
                    original_id=segment.original_id,
                    original_start=original_start,
                    original_end=original_end,
                    strand=strand,
                    reference_id=reference_id,
                    reference_start=qstart_i,
                    reference_end=qend_i,
                    reference_length=qend_i,
                    identity=float(attrs.get("Identity", "0")),
                    positive=float(attrs.get("Positive", "0")),
                    frameshifts=int(attrs.get("Frameshift", "0")),
                    inframe_stops=int(attrs.get("StopCodon", "0")),
                    rank=int(attrs.get("Rank", "999999")),
                    shard=segment.shard,
                    edge_distance=min(feature.start - 1, segment.length - feature.end),
                    features=[feature],
                    reference_hits={reference_id},
                )
            elif feature_type in {"CDS", "stop_codon"}:
                parent = attrs.get("Parent")
                if parent and parent in candidates:
                    candidates[parent].features.append(feature)
    return list(candidates.values())


def attach_reference_lengths(candidates: Iterable[Candidate], lengths: dict[str, int]) -> None:
    for candidate in candidates:
        candidate.reference_length = lengths.get(candidate.reference_id, candidate.reference_end)


def intervals_duplicate(left: Candidate, right: Candidate, threshold: float = 0.70) -> bool:
    if left.original_id != right.original_id or left.strand != right.strand:
        return False
    overlap = min(left.original_end, right.original_end) - max(left.original_start, right.original_start) + 1
    if overlap <= 0:
        return False
    return overlap / min(left.span, right.span) >= threshold


def deduplicate_candidates(candidates: list[Candidate]) -> list[Candidate]:
    # Build clusters around the best-scoring model at each locus.  A coordinate-
    # ordered connected-component pass is unsafe for tandem gene arrays: one long,
    # low-quality alignment can overlap two genuine neighbouring genes and merge
    # them transitively.  Quality-first anchors let that bridge contribute reference
    # evidence to one locus without ever merging two already distinct anchors.
    ordered = sorted(candidates, key=lambda item: item.quality_key(), reverse=True)
    clusters: list[list[Candidate]] = []
    for candidate in ordered:
        matches: list[tuple[float, float, int]] = []
        for index, cluster in enumerate(clusters):
            anchor = cluster[0]
            if not intervals_duplicate(candidate, anchor):
                continue
            overlap = (
                min(candidate.original_end, anchor.original_end)
                - max(candidate.original_start, anchor.original_start)
                + 1
            )
            reciprocal_overlap = overlap / max(candidate.span, anchor.span)
            shorter_overlap = overlap / min(candidate.span, anchor.span)
            matches.append((reciprocal_overlap, shorter_overlap, index))
        if matches:
            _, _, best_index = max(matches)
            clusters[best_index].append(candidate)
        else:
            clusters.append([candidate])
    kept: list[Candidate] = []
    for cluster in clusters:
        best = cluster[0]
        best.reference_hits = {
            reference_id for item in cluster for reference_id in item.reference_hits
        }
        kept.append(best)
    return sorted(kept, key=lambda item: (item.original_id, item.original_start, item.strand))


def extract_sequences(candidate: Candidate, virtual_sequence: str) -> None:
    cds_features = [feature for feature in candidate.features if feature.feature_type == "CDS"]
    cds_features.sort(key=lambda feature: feature.start, reverse=candidate.strand == "-")
    cds_parts: list[str] = []
    for feature in cds_features:
        part = virtual_sequence[feature.start - 1 : feature.end]
        if candidate.strand == "-":
            part = reverse_complement(part)
        cds_parts.append(part)
    # Internal exon phases describe how a codon is split across the splice
    # junction; deleting phase bases from every exon would corrupt such codons.
    # Only the leading phase of a partial 5' CDS is trimmed before translation.
    leading_phase = (
        int(cds_features[0].phase)
        if cds_features and cds_features[0].phase in {"0", "1", "2"}
        else 0
    )
    candidate.cds_sequence = "".join(cds_parts)[leading_phase:]
    span_start = min(feature.start for feature in candidate.features) - 1
    span_end = max(feature.end for feature in candidate.features)
    locus = virtual_sequence[span_start:span_end]
    candidate.locus_sequence = reverse_complement(locus) if candidate.strand == "-" else locus
    translated = translate(candidate.cds_sequence)
    candidate.has_start = candidate.cds_sequence.startswith("ATG")
    stop_features = [feature for feature in candidate.features if feature.feature_type == "stop_codon"]
    candidate.has_stop = bool(stop_features) or translated.endswith("*")
    candidate.internal_stops = translated.rstrip("*").count("*")
    candidate.protein_sequence = translated.rstrip("*")
    candidate.tm_helices = predict_tm_helices(candidate.protein_sequence)


HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3,
    "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5,
    "K": -3.9, "R": -4.5, "X": 0.0, "*": 0.0,
}


def predict_tm_helices(protein: str, window: int = 19, cutoff: float = 1.55) -> int:
    """Lightweight hydropathy evidence; this is deliberately not a TMHMM replacement."""
    if len(protein) < window:
        return 0
    hydrophobic_windows: list[int] = []
    for start in range(len(protein) - window + 1):
        score = sum(HYDROPATHY.get(amino_acid, 0.0) for amino_acid in protein[start : start + window]) / window
        charged = sum(amino_acid in "DEKR" for amino_acid in protein[start : start + window])
        if score >= cutoff and charged <= 2:
            hydrophobic_windows.append(start)
    if not hydrophobic_windows:
        return 0
    regions = 1
    previous = hydrophobic_windows[0]
    for start in hydrophobic_windows[1:]:
        if start - previous > window // 2:
            regions += 1
        previous = start
    return regions


def classify_candidate(
    candidate: Candidate,
    min_identity: float,
    min_report_coverage: float,
    min_protein_length: int,
    max_protein_length: int,
) -> None:
    length = len(candidate.protein_sequence)
    disrupted = candidate.frameshifts > 0 or candidate.inframe_stops > 0 or candidate.internal_stops > 0
    if candidate.coverage < min_report_coverage or candidate.identity < min_identity:
        candidate.status = "low_quality"
    elif disrupted:
        candidate.status = "pseudogene_candidate"
    elif length < min_protein_length:
        candidate.status = "partial"
    elif length > max_protein_length:
        candidate.status = "atypical_length"
    elif candidate.has_start and candidate.has_stop and candidate.coverage >= 0.70:
        candidate.status = "functional_candidate"
    else:
        candidate.status = "partial"

    length_ok = min_protein_length <= length <= max_protein_length
    tm_ok = 5 <= candidate.tm_helices <= 9
    if candidate.coverage >= 0.70 and candidate.identity >= 0.25 and length_ok and tm_ok:
        candidate.confidence = "high"
    elif candidate.coverage >= 0.40 and candidate.identity >= 0.18 and length >= 150 and 3 <= candidate.tm_helices <= 10:
        candidate.confidence = "medium"
    else:
        candidate.confidence = "low"

