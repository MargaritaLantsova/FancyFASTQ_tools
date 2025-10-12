#!/usr/bin/env python3
"""
Main script
"""

import argparse
from modules import sequence_tools, fastq_tools
from modules.fastq_tools import filter_fastq_stream

# Constants from module fastq_tools
DEFAULT_GC_BOUNDS = fastq_tools.DEFAULT_GC_BOUNDS
DEFAULT_LENGTH_BOUNDS = fastq_tools.DEFAULT_LENGTH_BOUNDS
DEFAULT_QUALITY_THRESHOLD = fastq_tools.DEFAULT_QUALITY_THRESHOLD


def run_dna_rna_tools(*args):
    """
    Function for DNA/RNA sequence operations.
    Supports both interactive mode and function call with arguments.
    """
    if args:
        return _process_sequences(*args)
    return _interactive_mode()


def _process_sequences(*args):
    """Process sequences and operation passed as arguments."""
    if len(args) < 2:
        raise ValueError("Sequence and operation must be provided")

    *sequences, operation = args

    for seq in sequences:
        if not seq:
            raise ValueError("Empty sequence provided")
        if not sequence_tools.is_valid_nucleic_acid(seq):
            raise ValueError(f"Invalid nucleic acid sequence: {seq}")

    if operation not in sequence_tools.OPERATIONS:
        raise ValueError(f"Unknown operation: {operation}")

    operation_func = sequence_tools.OPERATIONS[operation]
    results = [operation_func(seq) for seq in sequences]
    return results[0] if len(results) == 1 else results


def _interactive_mode():
    """Interactive mode for sequence operations."""
    print("=== DNA/RNA Tools ===")

    while True:
        print("\nChoose operation:")
        print("1 - Reverse sequence")
        print("2 - Complementary sequence")
        print("3 - Reverse complement sequence")
        print("4 - Transcribe DNA ↔ RNA")
        print("5 - Check sequence validity")
        print("6 - Determine type (DNA/RNA)")
        print("7 - GC content")
        print("0 - Exit")

        choice = input("Your choice: ").strip()

        if choice == "0":
            break
        if choice in ["1", "2", "3", "4", "5", "6", "7"]:
            sequence = input("Enter sequence: ").strip()
            if not sequence:
                print("Error: empty sequence entered")
                continue
            if choice != "5" and not sequence_tools.is_valid_nucleic_acid(sequence):
                print("Error: invalid nucleic acid sequence")
                continue

            try:
                operation_map = {
                    "1": ("reverse", sequence_tools.reverse),
                    "2": ("complement", sequence_tools.complement),
                    "3": (
                        "reverse complement",
                        sequence_tools.reverse_complement,
                    ),
                    "4": ("transcription", sequence_tools.transcribe),
                    "5": (
                        "validity check",
                        sequence_tools.is_valid_nucleic_acid,
                    ),
                    "6": (
                        "type determination",
                        lambda seq: (
                            "DNA"
                            if sequence_tools.is_dna(seq)
                            else (
                                "RNA"
                                if sequence_tools.is_rna(seq)
                                else "Unknown"
                            )
                        ),
                    ),
                    "7": ("GC content", sequence_tools.gc_content),
                }

                op_name, op_func = operation_map[choice]
                result = op_func(sequence)

                if choice == "7":
                    print(f"GC content: {result:.2f}%")
                elif choice == "6":
                    print(f"Type: {result}")
                elif choice == "5":
                    print(f"Valid sequence: {result}")
                else:
                    print(f"{op_name.capitalize()}: {result}")

            except ValueError as err:
                print(f"Error: {err}")
        else:
            print("Invalid choice, please try again")


def filter_fastq(
    seqs,
    gc_bounds=DEFAULT_GC_BOUNDS,
    length_bounds=DEFAULT_LENGTH_BOUNDS,
    quality_threshold=DEFAULT_QUALITY_THRESHOLD,
):
    """
    Filter FASTQ sequences by GC content, length and quality.
    """
    filtered_seqs = {}
    for seq_name, seq_data in seqs.items():
        passed = fastq_tools.filter_sequence(
            seq_data,
            gc_bounds,
            length_bounds,
            quality_threshold,
        )
        if passed:
            filtered_seqs[seq_name] = seq_data
    return filtered_seqs


def _cli():
    """Command-line interface."""
    parser = argparse.ArgumentParser(description="FancyFASTQ tools")
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
        help="Length bounds: one upper threshold or two values (min max)",
    )
    fastq_parser.add_argument(
        "--min-qual",
        type=int,
        default=0,
        help="Minimum average Phred+33",
    )

    args = parser.parse_args()

    if args.cmd == "fastq-filter":
        gc_bounds = (
            args.gc_bounds
            if len(args.gc_bounds) != 1
            else args.gc_bounds[0]
        )
        length_bounds = (
            args.length_bounds
            if len(args.length_bounds) != 1
            else args.length_bounds[0]
        )
        out_path, total, kept = filter_fastq_stream(
            input_fastq=args.input_fastq,
            output_fastq=args.output_fastq,
            gc_bounds=gc_bounds,
            length_bounds=length_bounds,
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


if __name__ == "__main__":
    _cli()
