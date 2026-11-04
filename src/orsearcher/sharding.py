from __future__ import annotations

import csv
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator, TextIO

from .fasta import open_text, write_record
from .utils import atomic_write_json, fingerprint, load_json


@dataclass
class Segment:
    virtual_id: str
    original_id: str
    original_description: str
    original_start: int  # 0-based, inclusive
    original_end: int  # 0-based, exclusive
    original_length: int
    shard: str
    split: bool

    @property
    def length(self) -> int:
        return self.original_end - self.original_start


@dataclass(frozen=True)
class ShardResult:
    directory: Path
    shards: tuple[Path, ...]
    mapping_path: Path
    segments: tuple[Segment, ...]
    total_bases: int
    sequence_count: int


class ShardWriter:
    def __init__(self, directory: Path, shard_size: int) -> None:
        self.directory = directory
        self.shard_size = shard_size
        self.index = 0
        self.current_bases = 0
        self.current_records = 0
        self.handle: TextIO | None = None
        self.path: Path | None = None
        self.paths: list[Path] = []

    def add(self, record_id: str, sequence: str) -> str:
        if self.handle is None or (self.current_records and self.current_bases + len(sequence) > self.shard_size):
            self._open_next()
        assert self.handle is not None and self.path is not None
        write_record(self.handle, record_id, sequence)
        self.current_bases += len(sequence)
        self.current_records += 1
        return self.path.name

    def _open_next(self) -> None:
        self.close_handle()
        self.index += 1
        self.path = self.directory / f"shard_{self.index:05d}.fa"
        self.handle = self.path.open("w", encoding="utf-8", newline="\n")
        self.paths.append(self.path)
        self.current_bases = 0
        self.current_records = 0

    def close_handle(self) -> None:
        if self.handle is not None:
            self.handle.close()
            self.handle = None

    def close(self) -> None:
        self.close_handle()


def stream_windows(genome: Path, chunk_size: int, overlap: int) -> Iterator[tuple[Segment, str]]:
    """Stream FASTA windows without loading a chromosome-sized record in memory."""
    step = chunk_size - overlap
    with open_text(genome) as handle:
        header: str | None = None
        record_id = ""
        buffer = bytearray()
        buffer_start = 0
        total_length = 0
        record_segments: list[Segment] = []
        seen_ids: set[str] = set()

        def finish_record() -> Iterator[tuple[Segment, str]]:
            nonlocal buffer
            if header is None:
                return
            if not record_segments:
                if not buffer:
                    raise ValueError(f"{genome}: record {record_id!r} has an empty sequence")
                segment = Segment("", record_id, header, 0, total_length, total_length, "", False)
                record_segments.append(segment)
                yield segment, buffer.decode("ascii")
            elif total_length > record_segments[-1].original_end:
                segment = Segment(
                    "", record_id, header, buffer_start, total_length, total_length, "", True
                )
                record_segments.append(segment)
                yield segment, buffer.decode("ascii")
            is_split = len(record_segments) > 1
            for segment in record_segments:
                segment.original_length = total_length
                segment.split = is_split
            buffer = bytearray()

        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                yield from finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{genome}:{line_no}: empty FASTA header")
                record_id = header.split()[0]
                if record_id in seen_ids:
                    raise ValueError(f"{genome}: duplicate FASTA ID {record_id!r}")
                seen_ids.add(record_id)
                buffer = bytearray()
                buffer_start = 0
                total_length = 0
                record_segments = []
                continue
            if header is None:
                raise ValueError(f"{genome}:{line_no}: sequence appears before FASTA header")
            sequence_line = "".join(line.split()).upper()
            try:
                buffer.extend(sequence_line.encode("ascii"))
            except UnicodeEncodeError as error:
                raise ValueError(f"{genome}:{line_no}: non-ASCII character in sequence") from error
            total_length += len(sequence_line)
            # Strictly greater preserves an exactly chunk-sized record as unsplit.
            while len(buffer) > chunk_size:
                segment = Segment(
                    "",
                    record_id,
                    header,
                    buffer_start,
                    buffer_start + chunk_size,
                    -1,
                    "",
                    True,
                )
                record_segments.append(segment)
                yield segment, buffer[:chunk_size].decode("ascii")
                buffer = buffer[step:]
                buffer_start += step
        if header is None:
            raise ValueError(f"{genome}: no FASTA records found")
        yield from finish_record()


def create_shards(
    genome: Path,
    directory: Path,
    chunk_size: int,
    overlap: int,
    shard_size: int,
    resume: bool = True,
) -> ShardResult:
    if not 0 <= overlap < chunk_size:
        raise ValueError("overlap must be >= 0 and smaller than chunk size")
    if shard_size < chunk_size:
        raise ValueError("shard size must be at least as large as chunk size")
    directory.mkdir(parents=True, exist_ok=True)
    metadata_path = directory / "shards.json"
    mapping_path = directory / "mapping.tsv"
    params = {"chunk_size": chunk_size, "overlap": overlap, "shard_size": shard_size}
    current_fingerprint = fingerprint([genome], params)
    if resume and metadata_path.exists() and mapping_path.exists():
        metadata = load_json(metadata_path)
        shards = tuple(directory / name for name in metadata.get("shards", []))
        if metadata.get("fingerprint") == current_fingerprint and shards and all(p.exists() for p in shards):
            segments = tuple(read_mapping(mapping_path))
            return ShardResult(
                directory,
                shards,
                mapping_path,
                segments,
                int(metadata["total_bases"]),
                int(metadata["sequence_count"]),
            )

    # Stale shard FASTAs are scoped to this work directory and are safe to replace.
    for stale in directory.glob("shard_*.fa"):
        stale.unlink()

    writer = ShardWriter(directory, shard_size)
    segments: list[Segment] = []
    sequence_count = 0
    virtual_number = 0
    try:
        last_original_id: str | None = None
        for segment, sequence in stream_windows(genome, chunk_size, overlap):
            if segment.original_id != last_original_id:
                sequence_count += 1
                last_original_id = segment.original_id
            virtual_number += 1
            virtual_id = f"ors{virtual_number:010d}"
            shard_name = writer.add(virtual_id, sequence)
            segment.virtual_id = virtual_id
            segment.shard = shard_name
            segments.append(segment)
    finally:
        writer.close()

    if not writer.paths:
        raise ValueError(f"{genome}: no sequences were written")
    total_bases = sum(
        segment.original_length
        for index, segment in enumerate(segments)
        if index == 0 or segment.original_id != segments[index - 1].original_id
    )
    write_mapping(mapping_path, segments)
    atomic_write_json(
        metadata_path,
        {
            "format_version": 1,
            "fingerprint": current_fingerprint,
            "genome": str(genome.resolve()),
            "parameters": params,
            "shards": [path.name for path in writer.paths],
            "segments": len(segments),
            "sequence_count": sequence_count,
            "total_bases": total_bases,
        },
    )
    return ShardResult(
        directory,
        tuple(writer.paths),
        mapping_path,
        tuple(segments),
        total_bases,
        sequence_count,
    )


MAPPING_FIELDS = [
    "virtual_id",
    "original_id",
    "original_description",
    "original_start",
    "original_end",
    "original_length",
    "shard",
    "split",
]


def write_mapping(path: Path, segments: list[Segment]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAPPING_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for segment in segments:
            row = asdict(segment)
            row["split"] = "1" if segment.split else "0"
            writer.writerow(row)
    os.replace(temporary, path)


def read_mapping(path: Path) -> Iterator[Segment]:
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield Segment(
                virtual_id=row["virtual_id"],
                original_id=row["original_id"],
                original_description=row["original_description"],
                original_start=int(row["original_start"]),
                original_end=int(row["original_end"]),
                original_length=int(row["original_length"]),
                shard=row["shard"],
                split=row["split"] == "1",
            )

