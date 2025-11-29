#!/usr/bin/env python3
"""
Entry points for sequence utilities and a small CLI for streaming FASTQ
filtering (HW5). This file intentionally exposes only two public functions:
- run_dna_rna_tools
- filter_fastq
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Union

from modules import sequence_tools, fastq_tools
from modules.fastq_tools import filter_fastq_stream

# Re-exported defaults used by callers / CLI
DEFAULT_GC_BOUNDS = fastq_tools.DEFAULT_GC_BOUNDS
DEFAULT_LENGTH_BOUNDS = fastq_tools.DEFAULT_LENGTH_BOUNDS
DEFAULT_QUALITY_THRESHOLD = fastq_tools.DEFAULT_QUALITY_THRESHOLD


def run_dna_rna_tools(*args: str) -> Union[str, List[str], None]:
    """
    High-level sequence operations.

    If positional arguments are provided, the last arg is an operation key
    (see sequence_tools.OPERATIONS), and the preceding args are sequences.
    Returns a single result for one sequence, or a list for many.

    If called without arguments, delegates to the interactive console
    implemented in modules.sequence_tools (keeping main.py minimal).
    """
    if not args:
        return sequence_tools.interactive_cli()

    if len(args) < 2:
        raise ValueError("Provide sequences and an operation key")

    *sequences, operation = args

    for seq in sequences:
        if not seq:
            raise ValueError("Empty sequence provided")
        if not sequence_tools.is_valid_nucleic_acid(seq):
            raise ValueError(f"Invalid nucleic acid sequence: {seq}")

    if operation not in sequence_tools.OPERATIONS:
        raise ValueError(f"Unknown operation: {operation}")

    func = sequence_tools.OPERATIONS[operation]
    results = [func(seq) for seq in sequences]
    return results[0] if len(results) == 1 else results


def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[int, Tuple[int, int]] = DEFAULT_GC_BOUNDS,
    length_bounds: Union[int, Tuple[int, int]] = DEFAULT_LENGTH_BOUNDS,
    quality_threshold: int = DEFAULT_QUALITY_THRESHOLD,
) -> Dict[str, Tuple[str, str]]:
    """
    In-memory FASTQ-like filtering for dicts:
    seqs: {read_name: (sequence, quality)} -> filtered dict with same shape.
    """
    filtered: Dict[str, Tuple[str, str]] = {}
    for name, (sequence, quality) in seqs.items():
        passed = fastq_tools.filter_sequence(
            (sequence, quality),
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
            quality_threshold=quality_threshold,
        )
        if passed:
            filtered[name] = (sequence, quality)
    return filtered


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FancyFASTQ tools (HW5)")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    fastq_parser = subparsers.add_parser(
        "fastq-filter",
        help="On-the-fly FASTQ filtering",
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
        help="Minimum average Phred+33",
    )

    args = parser.parse_args()

    if args.cmd == "fastq-filter":
        gc_arg: Union[int, Tuple[int, int]] = (
            args.gc_bounds if len(args.gc_bounds) != 1 else args.gc_bounds[0]
        )
        len_arg: Union[int, Tuple[int, int]] = (
            args.length_bounds
            if len(args.length_bounds) != 1
            else args.length_bounds[0]
        )
        out_path, total, kept = filter_fastq_stream(
            input_fastq=args.input_fastq,
            output_fastq=args.output_fastq,
            gc_bounds=gc_arg,
            length_bounds=len_arg,
            quality_threshold=args.min_qual,
        )
        print(
            "Input: {inp}\nOutput: {out}\nTotal: {tot}\nKept: {kept}".format(
                inp=args.input_fastq,
                out=out_path,
                tot=total,
                kept=kept,
            )
        )
