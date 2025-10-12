from __future__ import annotations
import os
import io
import gzip
import tempfile
import itertools
from typing import Generator, Iterable, Tuple, Optional

FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

def _open_maybe_gzip(path: str, mode: str = "rt") -> io.TextIOBase:
    if path.endswith(".gz"):
        # encoding only for text mode
        return gzip.open(path, mode, encoding="utf-8")  # type: ignore[arg-type]
    return open(path, mode, encoding="utf-8")  # type: ignore[call-arg]

def iter_fastq(path: str) -> Generator[Tuple[str, str, str], None, None]:
    """
    Streaming FASTQ: reads four-line records one by one.
    name without '@', seq and qual without '\n'. Validates seq/qual lengths.
    """
    with _open_maybe_gzip(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError(f"Ожидался '@' в заголовке, строка: {header!r}")
            name = header[1:].strip().split()[0]

            seq = fh.readline()
            if not seq:
                raise ValueError("Неожиданный конец файла после заголовка")
            seq = seq.strip()

            plus = fh.readline()
            if not plus or not plus.startswith("+"):
                raise ValueError("Отсутствует '+'-строка записи FASTQ")

            qual = fh.readline()
            if not qual:
                raise ValueError("Неожиданный конец файла в строке качества")
            qual = qual.strip()

            if len(seq) != len(qual):
                raise ValueError(
                    f"Длины seq({len(seq)}) и qual({len(qual)}) не совпадают для рида {name}"
                )
            yield name, seq, qual

def phred33_avg(qual: str) -> float:
    # Mean quality score (Phred+33)
    return sum(ord(c) - 33 for c in qual) / max(1, len(qual))

def gc_percent(seq: str) -> float:
    s = seq.upper()
    gc = sum(1 for b in s if b in ("G", "C"))
    return 100.0 * gc / max(1, len(s))

def _normalize_bounds(bounds, default_low: int, default_high: int) -> Tuple[int, int]:
    if isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        lo, hi = int(bounds[0]), int(bounds[1])
    elif isinstance(bounds, int):
        lo, hi = default_low, int(bounds)
    else:
        lo, hi = default_low, default_high
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi

def safe_filtered_path(output_fastq: str) -> str:
    """
    Returns a safe path under the filtered/ directory, with a unique filename.
    """
    os.makedirs("filtered", exist_ok=True)
    base = os.path.basename(output_fastq)
    if not base:
        base = "filtered.fastq"
    out_path = os.path.join("filtered", base)

    root, ext = os.path.splitext(out_path)
    # handle double extension .fastq.gz properly
    if base.endswith(".fastq.gz"):
        root = out_path[:-9]
        ext = ".fastq.gz"
    elif base.endswith(".fq.gz"):
        root = out_path[:-5]
        ext = ".fq.gz"

    i = 1
    candidate = out_path
    while os.path.exists(candidate):
        candidate = f"{root}__{i}{ext}"
        i += 1
    return candidate

def _atomic_write_text(final_path: str) -> Tuple[tempfile._TemporaryFileWrapper, str]:
    """
    Creates a temporary file next to final_path and returns (handle, tmp_path).
    After writing, the caller does os.replace(tmp_path, final_path).
    """
    dir_ = os.path.dirname(os.path.abspath(final_path)) or "."
    os.makedirs(dir_, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile(
        mode="wt", encoding="utf-8", prefix=".tmp_", dir=dir_, delete=False
    )
    return tmp, tmp.name

def write_fastq_record(handle: io.TextIOBase, name: str, seq: str, qual: str) -> None:
    handle.write(f"@{name}\n{seq}\n+\n{qual}\n")

def filter_fastq_stream(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: int = 0,
) -> Tuple[str, int, int]:
    """
    Reads input_fastq record-by-record, applies filters, and immediately writes
    passing reads into a safe file inside the filtered/ folder.
    Returns (output_path, total, kept).
    """
    # normalize thresholds
    gc_lo, gc_hi = _normalize_bounds(gc_bounds, 0, 100)
    len_lo, len_hi = _normalize_bounds(length_bounds, 0, 2**32)

    out_path = safe_filtered_path(output_fastq)
    tmp_handle, tmp_path = _atomic_write_text(out_path)
    kept = 0
    total = 0
    try:
        with tmp_handle:
            for name, seq, qual in iter_fastq(input_fastq):
                total += 1
                if not (len_lo <= len(seq) <= len_hi):
                    continue
                if not (gc_lo <= gc_percent(seq) <= gc_hi):
                    continue
                if phred33_avg(qual) < quality_threshold:
                    continue
                write_fastq_record(tmp_handle, name, seq, qual)
                kept += 1
        os.replace(tmp_path, out_path)  # atomic replace
    except Exception:
        # on exception — do not leave temp files behind
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        finally:
            raise
    return out_path, total, kept
