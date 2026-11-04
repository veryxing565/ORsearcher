from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable


SIZE_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([KMGT]?)(?:I?B)?\s*$", re.I)


def parse_size(value: str) -> int:
    match = SIZE_RE.match(value)
    if not match:
        raise ValueError(f"invalid size {value!r}; examples: 500M, 2G, 250000")
    number, unit = match.groups()
    multiplier = {"": 1, "K": 1000, "M": 1000**2, "G": 1000**3, "T": 1000**4}[unit.upper()]
    return int(float(number) * multiplier)


def sanitize_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return clean or "sample"


def fingerprint(paths: Iterable[Path], parameters: dict[str, Any]) -> str:
    digest = hashlib.sha256()
    for path in sorted((p.resolve() for p in paths), key=str):
        stat = path.stat()
        digest.update(str(path).encode())
        digest.update(str(stat.st_size).encode())
        digest.update(str(stat.st_mtime_ns).encode())
    digest.update(json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode())
    return digest.hexdigest()


def atomic_write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def load_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def executable(name_or_path: str) -> str | None:
    candidate = Path(name_or_path)
    if candidate.parent != Path(".") or candidate.is_absolute():
        return str(candidate.resolve()) if candidate.is_file() else None
    return shutil.which(name_or_path)


def command_version(command: str) -> str:
    attempts = ([command, "--version"], [command, "-V"], [command, "-h"], [command])
    for args in attempts:
        try:
            result = subprocess.run(args, text=True, capture_output=True, timeout=10, check=False)
        except (OSError, subprocess.TimeoutExpired):
            continue
        text = (result.stdout + "\n" + result.stderr).strip()
        if result.returncode == 0 and text:
            return text.splitlines()[0][:300]
    return "version unknown"


def setup_logging(log_path: Path, verbose: bool = False) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("orsearcher")
    for handler in logger.handlers:
        handler.close()
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s", "%Y-%m-%d %H:%M:%S")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger

