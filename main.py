#!/usr/bin/env python3
"""
Contains:
1) OOP classes for biological sequences (DNA/RNA/Protein)
2) FASTQ filtering using Biopython (length, mean quality, GC) + CLI
3) File utilities (FASTA multiline->oneline, BLAST best hits, GBK neighbors) + CLI

This script is intended to be the single entry point for the repository.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple, Union, overload
import argparse
import gzip
import os
import re


Index = Union[int, slice]

GCBounds = Union[int, Tuple[int, int]]
LengthBounds = Union[int, Tuple[int, int]]


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
    Implements: complement, reverse, reverse_complement (+ gc_content).
    Polymorphism is achieved via _complement_table() and _alphabet().
    """

    def complement(self) -> "NucleicAcidSequence":
        self.check_alphabet()
        trans = str.maketrans(self._complement_table())
        return self._new(self._seq.translate(trans))  # type: ignore[return-value]

    def reverse(self) -> "NucleicAcidSequence":
        self.check_alphabet()
        return self._new(self._seq[::-1])  # type: ignore[return-value]

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()

    def gc_content(self) -> float:
        """
        GC content percentage (0..100).
        """
        self.check_alphabet()
        if not self._seq:
            return 0.0
        gc = sum(1 for b in self._seq if b in ("G", "C"))
        return 100.0 * gc / len(self._seq)

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
# Task 2: FASTQ filtering with Biopython
# =========================

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
    return open(p, mode, encoding="utf-8", errors="replace")


def _ensure_filtered_dir() -> Path:
    out_dir = Path("filtered")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


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
# Merged from modules/bio_files_processor.py
# =========================

def ensure_unique_output(
    requested: Optional[str],
    default_name: str,
) -> str:
    """
    Build a unique output path. If the suggested path exists, append a
    numeric suffix before the extension.
    """
    base = default_name if not requested else (os.path.basename(requested) or default_name)
    root_dir = os.path.dirname(requested) if requested else ""
    out_dir = root_dir or "."
    os.makedirs(out_dir, exist_ok=True)

    root, ext = os.path.splitext(base)
    if not ext:
        ext = ".txt"

    candidate = os.path.join(out_dir, f"{root}{ext}")
    index = 1
    while os.path.exists(candidate):
        candidate = os.path.join(out_dir, f"{root}__{index}{ext}")
        index += 1
    return candidate


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: Optional[str] = None,
) -> str:
    """
    Read a FASTA where sequences may span multiple lines and write an output
    FASTA with one sequence per record line (header + single sequence line).
    """
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(input_fasta)

    if output_fasta is None:
        base = os.path.basename(input_fasta)
        root = re.sub(r"\.fa(sta)?(\.gz)?$", "", base, flags=re.IGNORECASE)
        output_fasta = f"{root}.oneline.fasta"

    out_path = ensure_unique_output(output_fasta, "oneline.fasta")

    def flush_record(
        header: Optional[str],
        seq_chunks: List[str],
        out_stream,
    ) -> None:
        if header is None:
            return
        seq = "".join(seq_chunks).replace(" ", "").replace("\t", "")
        out_stream.write(f">{header}\n{seq}\n")

    with open(input_fasta, "rt", encoding="utf-8") as in_file, open(
        out_path, "wt", encoding="utf-8"
    ) as out_file:
        header = None
        seq_parts: List[str] = []
        for line in in_file:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush_record(header, seq_parts, out_file)
                header = line[1:].strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
        flush_record(header, seq_parts, out_file)

    return out_path


_BLAST_SEC_RE = re.compile(
    r"^Sequences producing significant alignments:",
    re.IGNORECASE,
)
_SPLIT_COLS = re.compile(r"\s{2,}")


def parse_blast_output(input_file: str, output_file: str) -> str:
    """
    Expect a classic BLAST text output. For each section with the table
    'Sequences producing significant alignments:' take the first row and
    keep the Description column (first column). Write unique descriptions,
    sorted alphabetically, one per line.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(input_file)

    descriptions: List[str] = []
    with open(input_file, "rt", encoding="utf-8", errors="ignore") as in_file:
        lines = in_file.readlines()

    line_idx = 0
    while line_idx < len(lines):
        line = lines[line_idx]
        if _BLAST_SEC_RE.search(line):
            line_idx += 1
            while line_idx < len(lines):
                stripped = lines[line_idx].strip()
                is_header = stripped.lower().startswith("description")
                if stripped and not is_header:
                    break
                line_idx += 1
            if (
                line_idx < len(lines)
                and lines[line_idx].strip()
                and not lines[line_idx].startswith(">")
            ):
                row_text = lines[line_idx].rstrip("\n").strip()
                row = _SPLIT_COLS.split(row_text)
                if row:
                    descriptions.append(row[0])
        line_idx += 1

    uniq_sorted = sorted(set(descriptions), key=lambda s: s.lower())
    out_path = ensure_unique_output(
        output_file,
        "blast_best_hits.txt",
    )
    with open(out_path, "wt", encoding="utf-8") as out_file:
        for desc in uniq_sorted:
            out_file.write(desc + "\n")
    return out_path


def _parse_gbk_cds(input_gbk: str) -> List[Dict[str, str]]:
    """
    Naive FEATURES → CDS parser extracting /gene, /locus_tag
    and /translation.
    Order is preserved as in file (used as linear genomic order).
    """
    cds_list: List[Dict[str, str]] = []
    in_features = False
    current: Optional[Dict[str, str]] = None
    in_translation = False
    trans_buffer: List[str] = []

    with open(input_gbk, "rt", encoding="utf-8", errors="ignore") as in_file:
        for raw_line in in_file:
            line = raw_line.rstrip("\n")

            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue

            if re.match(r"^\s+CDS\s", line):
                if current:
                    if in_translation:
                        current["translation"] = (
                            "".join(trans_buffer).replace(" ", "").replace("\t", "")
                        )
                        in_translation = False
                        trans_buffer = []
                    cds_list.append(current)
                current = {"gene": "", "locus_tag": "", "translation": ""}
                continue

            if current is None:
                continue

            match_gene = re.match(r'^\s+/gene="([^"]+)"', line)
            if match_gene:
                current["gene"] = match_gene.group(1).strip()
                continue

            match_tag = re.match(r'^\s+/locus_tag="([^"]+)"', line)
            if match_tag:
                current["locus_tag"] = match_tag.group(1).strip()
                continue

            if re.match(r'^\s+/translation="', line):
                in_translation = True
                if '"/translation="' in line:
                    part = line.split('"/translation="')[-1]
                else:
                    part = line.split('/translation="', 1)[1]
                if part.endswith('"'):
                    in_translation = False
                    part = part[:-1]
                    current["translation"] = part.replace(" ", "").replace("\t", "")
                else:
                    trans_buffer = [part]
                continue

            if in_translation:
                stripped = line.strip()
                if stripped.endswith('"'):
                    trans_buffer.append(stripped[:-1])
                    current["translation"] = (
                        "".join(trans_buffer).replace(" ", "").replace("\t", "")
                    )
                    in_translation = False
                    trans_buffer = []
                else:
                    trans_buffer.append(stripped)

    if current:
        if in_translation:
            current["translation"] = (
                "".join(trans_buffer).replace(" ", "").replace("\t", "")
            )
        cds_list.append(current)

    return cds_list


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: Union[str, Sequence[str]],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "neighbors.fasta",
) -> str:
    """
    For each target gene (match by /gene or /locus_tag), write translations
    of n_before and n_after neighbor CDS entries into a FASTA file. Target
    genes themselves are not included.
    """
    if not os.path.exists(input_gbk):
        raise FileNotFoundError(input_gbk)

    if isinstance(genes, str):
        names = [g.strip() for g in re.split(r"[,\s;]+", genes) if g.strip()]
    else:
        names = [str(g).strip() for g in genes if str(g).strip()]
    names_lower = {g.lower() for g in names}

    cds = _parse_gbk_cds(input_gbk)

    def _name_of(record: Dict[str, str]) -> str:
        return (record.get("gene") or record.get("locus_tag") or "").strip()

    indices_of_targets = [
        idx for idx, record in enumerate(cds)
        if _name_of(record).lower() in names_lower
    ]
    if not indices_of_targets:
        raise ValueError("None of the genes of interest were found in the GBK.")

    indices_to_take = set()
    for target_index in indices_of_targets:
        start_left = max(0, target_index - n_before)
        for left_index in range(start_left, target_index):
            indices_to_take.add(left_index)

        end_right = min(len(cds), target_index + 1 + n_after)
        for right_index in range(target_index + 1, end_right):
            indices_to_take.add(right_index)

    out_path = ensure_unique_output(
        output_fasta,
        "neighbors.fasta",
    )
    with open(out_path, "wt", encoding="utf-8") as out_file:
        for idx in sorted(indices_to_take):
            record = cds[idx]
            name = _name_of(record) or f"CDS_{idx}"
            seq = (record.get("translation") or "").replace(" ", "")
            if not seq:
                continue
            header = f"{name}|idx={idx}"
            out_file.write(f">{header}\n{seq}\n")

    return out_path


# =========================
# Unified CLI
# =========================

def main() -> None:
    parser = argparse.ArgumentParser(description="Bioinformatics utilities (HW16)")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # FASTQ filter
    fastq_parser = subparsers.add_parser(
        "fastq-filter",
        help="FASTQ filtering (Biopython): length, mean quality, GC",
    )
    fastq_parser.add_argument("--input-fastq", required=True, help="Path to input .fastq[.gz]")
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

    # FASTA convert
    convert_parser = subparsers.add_parser(
        "convert-fasta",
        help="FASTA multiline to one-line",
    )
    convert_parser.add_argument("--input-fasta", required=True)
    convert_parser.add_argument("--output-fasta")

    # BLAST parse
    blast_parser = subparsers.add_parser(
        "parse-blast",
        help="BLAST txt → list of best-hit descriptions (sorted)",
    )
    blast_parser.add_argument("--input-file", required=True)
    blast_parser.add_argument("--output-file", required=True)

    # GBK neighbors
    gbk_parser = subparsers.add_parser(
        "gbk-neighbors",
        help="Pick neighbor CDS translations and write to FASTA",
    )
    gbk_parser.add_argument("--input-gbk", required=True)
    gbk_parser.add_argument(
        "--genes",
        required=True,
        nargs="+",
        help="Comma/space string or multiple args",
    )
    gbk_parser.add_argument("--n-before", type=int, default=1)
    gbk_parser.add_argument("--n-after", type=int, default=1)
    gbk_parser.add_argument("--output-fasta", default="neighbors.fasta")

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

    elif args.cmd == "convert-fasta":
        out_file = convert_multiline_fasta_to_oneline(
            args.input_fasta,
            args.output_fasta,
        )
        print(out_file)

    elif args.cmd == "parse-blast":
        out_file = parse_blast_output(
            args.input_file,
            args.output_file,
        )
        print(out_file)

    elif args.cmd == "gbk-neighbors":
        genes_arg: Union[str, Sequence[str]] = args.genes if len(args.genes) > 1 else args.genes[0]
        out_file = select_genes_from_gbk_to_fasta(
            input_gbk=args.input_gbk,
            genes=genes_arg,
            n_before=args.n_before,
            n_after=args.n_after,
            output_fasta=args.output_fasta,
        )
        print(out_file)


if __name__ == "__main__":
    main()
