from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


DNA_ALPHABET = set("ACGTUNRYKMSWBDHV-.")


@dataclass(frozen=True)
class FastaRecord:
    id: str
    description: str
    sequence: str


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if mode.startswith('r'):
        with path.open('rb') as probe:
            compressed = probe.read(2) == b'\x1f\x8b'
    else:
        compressed = path.suffix.lower() == '.gz'
    if compressed:
        return gzip.open(path, mode, encoding="utf-8")  # type: ignore[return-value]
    return path.open(mode, encoding="utf-8")


def read_fasta(path: Path) -> Iterator[FastaRecord]:
    """Read FASTA one record at a time; only one sequence is retained in memory."""
    with open_text(path) as handle:
        header: str | None = None
        chunks: list[str] = []
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _make_record(header, chunks, path)
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_no}: sequence appears before FASTA header")
                chunks.append("".join(line.split()))
        if header is not None:
            yield _make_record(header, chunks, path)
        elif not chunks:
            raise ValueError(f"{path}: no FASTA records found")


def _make_record(header: str, chunks: list[str], path: Path) -> FastaRecord:
    if not header:
        raise ValueError(f"{path}: empty FASTA header")
    sequence = "".join(chunks).upper()
    if not sequence:
        raise ValueError(f"{path}: record {header!r} has an empty sequence")
    return FastaRecord(header.split()[0], header, sequence)


def write_record(handle: TextIO, record_id: str, sequence: str, description: str = "") -> None:
    header = record_id if not description else f"{record_id} {description}"
    handle.write(f">{header}\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset : offset + 80] + "\n")


def write_fasta(path: Path, records: Iterable[FastaRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            write_record(handle, record.id, record.sequence, record.description)


def detect_sequence_type(sequence: str) -> str:
    cleaned = sequence.upper().replace("*", "")
    if not cleaned:
        raise ValueError("empty sequence")
    dna_fraction = sum(base in DNA_ALPHABET for base in cleaned) / len(cleaned)
    return "dna" if dna_fraction >= 0.90 else "protein"


COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


CODON_TABLE_1 = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(sequence: str, trim_partial: bool = True) -> str:
    sequence = sequence.upper().replace("U", "T")
    if trim_partial:
        sequence = sequence[: len(sequence) - len(sequence) % 3]
    return "".join(CODON_TABLE_1.get(sequence[i : i + 3], "X") for i in range(0, len(sequence) - 2, 3))


def longest_orf(sequence: str) -> tuple[str, str]:
    """Return (protein, strand/frame label) for the best six-frame ORF.

    Complete ATG-starting segments are preferred, but incomplete reference CDSs
    remain usable.
    """
    candidates: list[tuple[tuple[int, int], str, str]] = []
    for strand, dna in (("+", sequence), ("-", reverse_complement(sequence))):
        for frame in range(3):
            protein = translate(dna[frame:])
            for segment in protein.split("*"):
                if segment:
                    methionine = segment.find("M")
                    if methionine >= 0:
                        complete = segment[methionine:]
                        candidates.append(((1, len(complete)), complete, f"{strand}{frame + 1}"))
                    candidates.append(((0, len(segment)), segment, f"{strand}{frame + 1}"))
    if not candidates:
        return "", "+1"
    _, protein, label = max(candidates, key=lambda item: item[0])
    return protein, label

