#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from typing import Dict, List, Optional

# ---------- Safe-write utilities ------------------------------------

def ensure_unique_output(
    requested: Optional[str],
    default_name: str,
) -> str:
    """
    Build a unique output path. If the suggested path exists, append a
    numeric suffix before the extension.
    """
    base = default_name if not requested else (
        os.path.basename(requested) or default_name
    )
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


# ---------- 1) FASTA: multiline → one-line ------------------------------

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
                # new header — flush previous sequence
                flush_record(header, seq_parts, out_file)
                header = line[1:].strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
        # final flush
        flush_record(header, seq_parts, out_file)

    return out_path


# ---------- 2) BLAST: best hit per QUERY --------------------------------

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
            # skip possible headers / empty lines
            while line_idx < len(lines):
                stripped = lines[line_idx].strip()
                is_header = stripped.lower().startswith("description")
                if stripped and not is_header:
                    break
                line_idx += 1
            # first meaningful row is the best hit
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


# ---------- 3) GBK: neighbors of genes of interest to FASTA --------------

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

            # start of CDS
            if re.match(r"^\s+CDS\s", line):
                if current:
                    # finalize previous CDS (if translation still open)
                    if in_translation:
                        current["translation"] = (
                            "".join(trans_buffer)
                            .replace(" ", "")
                            .replace("\t", "")
                        )
                        in_translation = False
                        trans_buffer = []
                    cds_list.append(current)
                current = {"gene": "", "locus_tag": "", "translation": ""}
                continue

            if current is None:
                # ignore non-CDS features
                continue

            # qualifiers parsing
            match_gene = re.match(r'^\s+/gene="([^"]+)"', line)
            if match_gene:
                current["gene"] = match_gene.group(1).strip()
                continue

            match_tag = re.match(r'^\s+/locus_tag="([^"]+)"', line)
            if match_tag:
                current["locus_tag"] = match_tag.group(1).strip()
                continue

            # /translation may span multiple lines
            if re.match(r'^\s+/translation="', line):
                in_translation = True
                # capture everything after the opening quote
                if '"/translation="' in line:
                    part = line.split('"/translation="')[-1]
                else:
                    part = line.split('/translation="', 1)[1]
                # closing quote on the same line?
                if part.endswith('"'):
                    in_translation = False
                    part = part[:-1]
                    current["translation"] = (
                        part.replace(" ", "").replace("\t", "")
                    )
                else:
                    trans_buffer = [part]
                continue

            if in_translation:
                # look for closing quote
                stripped = line.strip()
                if stripped.endswith('"'):
                    trans_buffer.append(stripped[:-1])
                    current["translation"] = (
                        "".join(trans_buffer)
                        .replace(" ", "")
                        .replace("\t", "")
                    )
                    in_translation = False
                    trans_buffer = []
                else:
                    trans_buffer.append(stripped)

    # finalize the last CDS
    if current:
        if in_translation:
            current["translation"] = (
                "".join(trans_buffer).replace(" ", "").replace("\t", "")
            )
        cds_list.append(current)

    return cds_list


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes,
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
        raise ValueError(
            "None of the genes of interest were found in the GBK."
        )

    indices_to_take = set()
    for target_index in indices_of_targets:
        # left block
        start_left = max(0, target_index - n_before)
        for left_index in range(start_left, target_index):
            indices_to_take.add(left_index)
        # right block
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


# ---------------- CLI ----------------------------------------------

def _cli() -> None:
    """Command-line interface for the bio files processor."""
    parser = argparse.ArgumentParser(description="Bio files processor")
    subparsers = parser.add_subparsers(dest="cmd", required=True)

    convert_parser = subparsers.add_parser(
        "convert-fasta",
        help="FASTA multiline to one-line",
    )
    convert_parser.add_argument("--input-fasta", required=True)
    convert_parser.add_argument("--output-fasta")

    blast_parser = subparsers.add_parser(
        "parse-blast",
        help="BLAST txt → list of best-hit descriptions (sorted)",
    )
    blast_parser.add_argument("--input-file", required=True)
    blast_parser.add_argument("--output-file", required=True)

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

    if args.cmd == "convert-fasta":
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
        genes_arg = args.genes if len(args.genes) > 1 else args.genes[0]
        out_file = select_genes_from_gbk_to_fasta(
            input_gbk=args.input_gbk,
            genes=genes_arg,
            n_before=args.n_before,
            n_after=args.n_after,
            output_fasta=args.output_fasta,
        )
        print(out_file)


if __name__ == "__main__":
    _cli()
