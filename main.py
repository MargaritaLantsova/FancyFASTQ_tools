#!/usr/bin/env python3
"""
Contains:
- OOP classes for biological sequences (DNA/RNA/Protein)
- FASTQ filtering using Biopython (length, mean quality, GC)
- Small CLI command: fastq-filter
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Iterator, Optional, Sequence, Tuple, Union, overload
import gzip


Index = Union[int, slice]


# =========================
# Abstract sequences
# =========================

class BiologicalSequence(ABC):
    """
    Abstract interface for biological sequences.
    Required:
    - len(seq)
    - indexing and slicing
    - pretty print
    - alphabet validation
    """

    def __init__(self, sequence: str) -> None:
        if not isinstance(sequence, str):
            raise TypeError("sequence must be a str")
        if not sequence:
            raise ValueError("Empty sequence provided")
        self._seq: str = sequence.upper()

    def __len__(self) -> int:
        return len(self._seq)

    @overload
    def __getitem__(self, idx: int) -> str: ...
    @overload
    def __getitem__(self, idx: slice) -> "BiologicalSequence": ...

    def __getitem__(self, idx: Index) -> Union[str, "BiologicalSequence"]:
        if isinstance(idx, int):
            return self._seq[idx]
        if isinstance(idx, slice):
            return self._new(self._seq[idx])
        raise TypeError("Index must be int or slice")

    def __iter__(self) -> Iterator[str]:
        return iter(self._seq)

    def __str__(self) -> str:
        preview = self._seq if len(self._seq) <= 60 else (self._seq[:57] + "...")
        return f"{self.__class__.__name__}(len={len(self)}): {preview}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._seq!r})"

    @property
    def sequence(self) -> str:
        return self._seq

    def check_alphabet(self) -> None:
        allowed = self._alphabet()
        bad = sorted({ch for ch in self._seq if ch not in allowed})
        if bad:
            raise ValueError(
                f"{self.__class__.__name__}: invalid symbols {bad}. Allowed: {sorted(allowed)}"
            )

    @abstractmethod
    def _alphabet(self) -> set[str]:
        raise NotImplementedError

    @abstractmethod
    def _new(self, seq: str) -> "BiologicalSequence":
        raise NotImplementedError


class NucleicAcidSequence(BiologicalSequence, ABC):
    """
    Parent for DNA/RNA.
    Implements: complement, reverse, reverse_complement
    Polymorphism is achieved via _complement_table() and _alphabet().
    """
    def gc_content(self) -> float:
        """
        GC content percentage (0..100).
        """
        self.check_alphabet()
        if not self._seq:
            return 0.0
        gc = sum(1 for b in self._seq if b in ("G", "C"))
        return 100.0 * gc / len(self._seq)
        
    def complement(self) -> "NucleicAcidSequence":
        self.check_alphabet()
        trans = str.maketrans(self._complement_table())
        return self._new(self._seq.translate(trans))  # type: ignore[return-value]

    def reverse(self) -> "NucleicAcidSequence":
        self.check_alphabet()
        return self._new(self._seq[::-1])  # type: ignore[return-value]

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()

    @abstractmethod
    def _complement_table(self) -> Dict[str, str]:
        raise NotImplementedError


class DNASequence(NucleicAcidSequence):
    def _alphabet(self) -> set[str]:
        return {"A", "T", "G", "C", "N"}

    def _complement_table(self) -> Dict[str, str]:
        return {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    def _new(self, seq: str) -> "DNASequence":
        return DNASequence(seq)

    def transcribe(self) -> "RNASequence":
        self.check_alphabet()
        return RNASequence(self._seq.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    def _alphabet(self) -> set[str]:
        return {"A", "U", "G", "C", "N"}

    def _complement_table(self) -> Dict[str, str]:
        return {"A": "U", "U": "A", "G": "C", "C": "G", "N": "N"}

    def _new(self, seq: str) -> "RNASequence":
        return RNASequence(seq)


class AminoAcidSequence(BiologicalSequence):
    """
    Protein sequence.
    Includes one domain-sensible method: approximate molecular weight.
    """

    def _alphabet(self) -> set[str]:
        # 20 aa + common ambiguous/rare + stop
        return set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")

    def _new(self, seq: str) -> "AminoAcidSequence":
        return AminoAcidSequence(seq)

    def molecular_weight_approx(self) -> float:
        self.check_alphabet()
        masses = {
            "A": 71.0788, "C": 103.1388, "D": 115.0886, "E": 129.1155, "F": 147.1766,
            "G": 57.0519, "H": 137.1411, "I": 113.1594, "K": 128.1741, "L": 113.1594,
            "M": 131.1926, "N": 114.1038, "P": 97.1167, "Q": 128.1307, "R": 156.1875,
            "S": 87.0782, "T": 101.1051, "V": 99.1326, "W": 186.2132, "Y": 163.1760,
        }
        total = 0.0
        for aa in self._seq:
            if aa in masses:
                total += masses[aa]
            else:
                raise ValueError(
                    f"Cannot estimate MW: unsupported/ambiguous amino acid '{aa}'. "
                    "Remove/resolve ambiguous symbols (X/B/Z/J/U/O/*) first."
                )
        return total + 18.01528


# =========================
# FASTQ filtering with Biopython
# =========================

GCBounds = Union[int, Tuple[int, int]]
LengthBounds = Union[int, Tuple[int, int]]

DEFAULT_GC_BOUNDS: GCBounds = (0, 100)
DEFAULT_LENGTH_BOUNDS: LengthBounds = (0, 2**32)
DEFAULT_QUALITY_THRESHOLD: int = 0


def _normalize_bounds(bounds: Union[int, Tuple[int, int], Sequence[int]]) -> Tuple[int, int]:
    if isinstance(bounds, int):
        return (0, bounds)
    if isinstance(bounds, tuple) and len(bounds) == 2:
        return (int(bounds[0]), int(bounds[1]))
    if isinstance(bounds, Sequence):
        b = list(bounds)
        if len(b) == 1:
            return (0, int(b[0]))
        if len(b) == 2:
            return (int(b[0]), int(b[1]))
    raise ValueError("Bounds must be int, (min,max), or a sequence of 1-2 ints")


def _open_maybe_gz(path: Union[str, Path], mode: str):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, mode)
    return open(p, mode, encoding="utf-8", errors="replace")  # text mode


def filter_fastq(
    input_fastq: Union[str, Path],
    output_fastq: Union[str, Path],
    gc_bounds: GCBounds = DEFAULT_GC_BOUNDS,
    length_bounds: LengthBounds = DEFAULT_LENGTH_BOUNDS,
    quality_threshold: int = DEFAULT_QUALITY_THRESHOLD,
) -> Tuple[str, int, int]:
    """
    FASTQ filtering using Biopython.

    Filters reads by:
    - length_bounds: int (0..max) or (min,max)
    - quality_threshold: minimum mean Phred quality (integer)
    - gc_bounds: int (0..max) or (min,max), in percent

    Returns: (output_path, total_reads, kept_reads)
    """
    try:
        from Bio import SeqIO
        from Bio.SeqUtils import gc_fraction
    except ImportError as e:
        raise ImportError("Biopython is required. Add 'biopython' to requirements.txt") from e

    gc_min, gc_max = _normalize_bounds(gc_bounds)
    len_min, len_max = _normalize_bounds(length_bounds)

    in_path = Path(input_fastq)
    out_path = Path(output_fastq)

    total = 0
    kept = 0

    # SeqIO.parse wants a handle; for gz use gzip in text mode ("rt")
    with _open_maybe_gz(in_path, "rt") as hin, _open_maybe_gz(out_path, "wt") as hout:
        out_records = []
        for rec in SeqIO.parse(hin, "fastq"):
            total += 1

            seq_len = len(rec.seq)
            if seq_len < len_min or seq_len > len_max:
                continue

            quals = rec.letter_annotations.get("phred_quality")
            if not quals:
                continue
            mean_q = sum(quals) / len(quals)
            if mean_q < quality_threshold:
                continue

            gc_percent = gc_fraction(str(rec.seq).upper()) * 100.0
            if gc_percent < gc_min or gc_percent > gc_max:
                continue

            out_records.append(rec)
            kept += 1

        SeqIO.write(out_records, hout, "fastq")

    return (str(out_path), total, kept)


# =========================
# CLI
# =========================

def _ensure_filtered_dir() -> Path:
    out_dir = Path("filtered")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FancyFASTQ tools (HW16)")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    fastq_parser = subparsers.add_parser(
        "fastq-filter",
        help="On-the-fly FASTQ filtering (Biopython)",
    )
    fastq_parser.add_argument(
        "--input-fastq",
        required=True,
        help="Path to input .fastq[.gz]",
    )
    fastq_parser.add_argument(
        "--output-fastq",
        required=True,
        help="Output FASTQ name (written to ./filtered)",
    )
    fastq_parser.add_argument(
        "--gc-bounds",
        nargs="+",
        type=int,
        default=[0, 100],
        help="GC bounds: one upper threshold or two values (min max)",
    )
    fastq_parser.add_argument(
        "--length-bounds",
        nargs="+",
        type=int,
        default=[0, 2**32],
        help="Length bounds: one upper or two values (min max)",
    )
    fastq_parser.add_argument(
        "--min-qual",
        type=int,
        default=0,
        help="Minimum mean Phred quality (integer)",
    )

    args = parser.parse_args()

    if args.cmd == "fastq-filter":
        gc_arg: GCBounds = args.gc_bounds if len(args.gc_bounds) != 1 else args.gc_bounds[0]
        len_arg: LengthBounds = (
            args.length_bounds if len(args.length_bounds) != 1 else args.length_bounds[0]
        )

        out_dir = _ensure_filtered_dir()
        out_path = out_dir / args.output_fastq

        out_path_str, total, kept = filter_fastq(
            input_fastq=args.input_fastq,
            output_fastq=out_path,
            gc_bounds=gc_arg,
            length_bounds=len_arg,
            quality_threshold=args.min_qual,
        )

        print(
            "Input: {inp}\nOutput: {out}\nTotal: {tot}\nKept: {kept}".format(
                inp=args.input_fastq,
                out=out_path_str,
                tot=total,
                kept=kept,
            )
        )
