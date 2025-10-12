from __future__ import annotations
import os
import io
import gzip
import tempfile
from typing import Generator, Tuple

FASTQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def _open_maybe_gzip(path: str, mode: str = "rt") -> io.TextIOBase:
    if path.endswith(".gz"):
        # encoding only in text mode
        return gzip.open(path, mode, encoding="utf-8")  # type: ignore[arg-type]
    return open(path, mode, encoding="utf-8")  # type: ignore[call-arg]


def iter_fastq(path: str) -> Generator[Tuple[str, str, str], None, None]:
    """
    Yield FASTQ records one-by-one as (name, seq, qual).
    `name` is without '@'; `seq` and `qual` are newline-stripped.
    """
    with _open_maybe_gzip(path, "rt") as file_in:
        while True:
            header = file_in.readline()
            if not header:
                break
            if not header.startswith("@"):
                raise ValueError(
                    "Expected '@' at record header, got: {!r}".format(header)
                )
            name = header[1:].strip().split()[0]

            seq_line = file_in.readline()
            if not seq_line:
                raise ValueError("Unexpected EOF after header")
            seq = seq_line.strip()

            plus = file_in.readline()
            if not plus or not plus.startswith("+"):
                raise ValueError("Missing '+' line in FASTQ record")

            qual_line = file_in.readline()
            if not qual_line:
                raise ValueError("Unexpected EOF in quality line")
            qual = qual_line.strip()

            if len(seq) != len(qual):
                raise ValueError(
                    "Length mismatch: seq({}) vs qual({}) for read {}".format(
                        len(seq), len(qual), name
                    )
                )
            yield name, seq, qual


def phred33_avg(qual: str) -> float:
    """Mean quality score (Phred+33)."""
    return sum(ord(ch) - 33 for ch in qual) / max(1, len(qual))


def gc_percent(seq: str) -> float:
    """GC percentage in the sequence."""
    seq_upper = seq.upper()
    gc_count = sum(1 for base in seq_upper if base in ("G", "C"))
    return 100.0 * gc_count / max(1, len(seq_upper))


def normalize_bounds_pair(
    bounds, default_low: int, default_high: int
) -> Tuple[int, int]:
    """Normalize bounds to a (low, high) integer pair."""
    if isinstance(bounds, (tuple, list)) and len(bounds) == 2:
        low, high = int(bounds[0]), int(bounds[1])
    elif isinstance(bounds, int):
        low, high = default_low, int(bounds)
    else:
        low, high = default_low, default_high
    if low > high:
        low, high = high, low
    return low, high


def safe_filtered_path(output_fastq: str) -> str:
    """
    Return a unique path under ./filtered for the requested output filename.
    """
    os.makedirs("filtered", exist_ok=True)
    base = os.path.basename(output_fastq) or "filtered.fastq"
    out_path = os.path.join("filtered", base)

    root, ext = os.path.splitext(out_path)
    # handle double extensions properly
    if base.endswith(".fastq.gz"):
        root = out_path[:-9]
        ext = ".fastq.gz"
    elif base.endswith(".fq.gz"):
        root = out_path[:-5]
        ext = ".fq.gz"

    index = 1
    candidate = out_path
    while os.path.exists(candidate):
        candidate = f"{root}__{index}{ext}"
        index += 1
    return candidate


def _atomic_write_text(final_path: str) -> Tuple[tempfile._TemporaryFileWrapper, str]:
    """
    Create a temp file next to `final_path` and return (handle, tmp_path).
    After writing, replace tmp with the final path atomically.
    """
    dir_path = os.path.dirname(os.path.abspath(final_path)) or "."
    os.makedirs(dir_path, exist_ok=True)
    tmp = tempfile.NamedTemporaryFile(
        mode="wt", encoding="utf-8", prefix=".tmp_", dir=dir_path, delete=False
    )
    return tmp, tmp.name


def write_fastq_record(
    out_stream: io.TextIOBase, name: str, seq: str, qual: str
) -> None:
    """Write a single FASTQ record to an open text stream."""
    out_stream.write(f"@{name}\n{seq}\n+\n{qual}\n")


def filter_fastq_stream(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold: int = 0,
) -> Tuple[str, int, int]:
    """
    Read records from `input_fastq`, filter on the fly, and write passing reads
    to a safe file under ./filtered. Return (output_path, total, kept).
    """
    gc_low, gc_high = normalize_bounds_pair(gc_bounds, 0, 100)
    len_low, len_high = normalize_bounds_pair(length_bounds, 0, 2**32)

    out_path = safe_filtered_path(output_fastq)
    tmp_handle, tmp_path = _atomic_write_text(out_path)
    kept = 0
    total = 0
    try:
        with tmp_handle:
            for name, seq, qual in iter_fastq(input_fastq):
                total += 1
                if not (len_low <= len(seq) <= len_high):
                    continue
                if not (gc_low <= gc_percent(seq) <= gc_high):
                    continue
                if phred33_avg(qual) < quality_threshold:
                    continue
                write_fastq_record(tmp_handle, name, seq, qual)
                kept += 1
        os.replace(tmp_path, out_path)  # atomic replace
    except Exception:
        # make sure no temp files are left behind
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        finally:
            raise
    return out_path, total, kept
