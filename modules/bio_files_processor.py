#!/usr/bin/env python3
from __future__ import annotations
import os
import re
import argparse
import itertools
from typing import List, Tuple, Dict, Optional

# ---------- Safe-write utilities ----------
def _safe_out_path(path: Optional[str], default_name: str) -> str:
    base = default_name if not path else os.path.basename(path) or default_name
    root_dir = os.path.dirname(path) if path else ""
    out_dir = root_dir if root_dir else "."
    os.makedirs(out_dir, exist_ok=True)
    root, ext = os.path.splitext(base)
    if not ext:
        ext = ".txt"
    candidate = os.path.join(out_dir, f"{root}{ext}")
    i = 1
    while os.path.exists(candidate):
        candidate = os.path.join(out_dir, f"{root}__{i}{ext}")
        i += 1
    return candidate

# ---------- 1) FASTA: multiline → one-line ----------
def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: Optional[str] = None) -> str:
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(input_fasta)
    # by default — next to the source, with .oneline.fasta suffix
    if output_fasta is None:
        base = os.path.basename(input_fasta)
        root = re.sub(r"\.fa(sta)?(\.gz)?$", "", base, flags=re.IGNORECASE)
        output_fasta = f"{root}.oneline.fasta"
    out_path = _safe_out_path(output_fasta, "oneline.fasta")

    def flush(header: Optional[str], seq_chunks: List[str], fh_out):
        if header is None:
            return
        seq = "".join(seq_chunks).replace(" ", "").replace("\t", "")
        fh_out.write(f">{header}\n{seq}\n")

    with open(input_fasta, "rt", encoding="utf-8") as fh_in, \
         open(out_path, "wt", encoding="utf-8") as fh_out:
        header = None
        seq_parts: List[str] = []
        for line in fh_in:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # new header — flush previous sequence
                flush(header, seq_parts, fh_out)
                header = line[1:].strip()
                seq_parts = []
            elif line.strip():
                seq_parts.append(line.strip())
        # final flush
        flush(header, seq_parts, fh_out)
    return out_path

# ---------- 2) BLAST: best hit per QUERY ----------
_BLAST_SEC_RE = re.compile(r"^Sequences producing significant alignments:", re.I)
_SPLIT_COLS = re.compile(r"\s{2,}")

def parse_blast_output(input_file: str, output_file: str) -> str:
    """
    Expect a classic BLAST plain-text output.
    For each section with the table "Sequences producing significant alignments:",
    take the FIRST non-empty row of that table and grab its Description (first column).
    Result: a sorted set of unique descriptions, one per line.
    """
    if not os.path.exists(input_file):
        raise FileNotFoundError(input_file)
    hits = []
    with open(input_file, "rt", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        line = lines[i]
        if _BLAST_SEC_RE.search(line):
            # usually the next line is empty or a header/separator
            i += 1
            # skip possible headers/empty lines
            while i < len(lines) and (not lines[i].strip() or lines[i].strip().lower().startswith("description")):
                i += 1
            # the first meaningful row of the table — that's the best hit
            if i < len(lines) and lines[i].strip() and not lines[i].startswith(">"):
                row = _SPLIT_COLS.split(lines[i].rstrip("\n").strip())
                if row:
                    hits.append(row[0])
        i += 1

    uniq_sorted = sorted(set(hits), key=lambda s: s.lower())
    out_path = _safe_out_path(output_file, "blast_best_hits.txt")
    with open(out_path, "wt", encoding="utf-8") as fh:
        for h in uniq_sorted:
            fh.write(h + "\n")
    return out_path

# ---------- 3) GBK: neighbors of genes of interest to FASTA ----------
def _parse_gbk_cds(input_gbk: str) -> List[Dict[str, str]]:
    """
    Naive FEATURES → CDS parser extracting /gene, /locus_tag, and /translation.
    Order is preserved as in file (used as a linear genomic order).
    """
    cds_list: List[Dict[str, str]] = []
    in_features = False
    current: Optional[Dict[str, str]] = None
    in_translation = False
    tr_buf: List[str] = []

    with open(input_gbk, "rt", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if line.startswith("FEATURES"):
                in_features = True
                continue
            if not in_features:
                continue
            # start of CDS
            if re.match(r"^\s+CDS\s", line):
                if current:
                    # finalize previous CDS (if translation was still open)
                    if in_translation:
                        current["translation"] = "".join(tr_buf).replace(" ", "").replace("\t", "")
                        in_translation = False
                        tr_buf = []
                    cds_list.append(current)
                current = {"gene": "", "locus_tag": "", "translation": ""}
                continue

            if current is None:
                # ignore non-CDS features
                continue

            # qualifiers parsing
            m = re.match(r'^\s+/gene="([^"]+)"', line)
            if m:
                current["gene"] = m.group(1).strip()
                continue
            m = re.match(r'^\s+/locus_tag="([^"]+)"', line)
            if m:
                current["locus_tag"] = m.group(1).strip()
                continue

            # /translation can span multiple lines
            if re.match(r'^\s+/translation="', line):
                in_translation = True
                # capture everything after the opening quote
                part = line.split('"/translation="')[-1] if '"/translation="' in line else line.split('/translation="',1)[1]
                # closing quote on the same line?
                if part.endswith('"'):
                    in_translation = False
                    part = part[:-1]
                    current["translation"] = part.replace(" ", "").replace("\t", "")
                else:
                    tr_buf = [part]
                continue

            if in_translation:
                # look for closing quote
                if line.strip().endswith('"'):
                    tr_buf.append(line.strip()[:-1])
                    current["translation"] = "".join(tr_buf).replace(" ", "").replace("\t", "")
                    in_translation = False
                    tr_buf = []
                else:
                    tr_buf.append(line.strip())

    # finalize the last CDS
    if current:
        if in_translation:
            current["translation"] = "".join(tr_buf).replace(" ", "").replace("\t", "")
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
    Take the CDS list in file order. For each gene of interest (match on /gene or /locus_tag),
    save n_before and n_after neighbors (DO NOT include the genes of interest themselves).
    """
    if not os.path.exists(input_gbk):
        raise FileNotFoundError(input_gbk)

    # normalize the list of names
    if isinstance(genes, str):
        # allow a single string like "acrA, tolC  ;  marA"
        names = [g.strip() for g in re.split(r"[,\s;]+", genes) if g.strip()]
    else:
        names = [str(g).strip() for g in genes if str(g).strip()]
    names_lower = {g.lower() for g in names}

    cds = _parse_gbk_cds(input_gbk)
    # name indexer
    def _name_of(rec: Dict[str, str]) -> str:
        return (rec.get("gene") or rec.get("locus_tag") or "").strip()

    indices_of_targets = [
        i for i, rec in enumerate(cds)
        if _name_of(rec).lower() in names_lower
    ]

    if not indices_of_targets:
        raise ValueError("None of the genes of interest were found in the GBK.")

    to_take = set()
    for idx in indices_of_targets:
        # left block
        for j in range(max(0, idx - n_before), idx):
            to_take.add(j)
        # right block
        for j in range(idx + 1, min(len(cds), idx + 1 + n_after)):
            to_take.add(j)

    # write FASTA (one-line sequences)
    out_path = _safe_out_path(output_fasta, "neighbors.fasta")
    with open(out_path, "wt", encoding="utf-8") as fh:
        for i in sorted(to_take):
            rec = cds[i]
            name = _name_of(rec) or f"CDS_{i}"
            seq = (rec.get("translation") or "").replace(" ", "")
            if not seq:
                # skip CDS without translation
                continue
            header = f"{name}|idx={i}"
            fh.write(f">{header}\n{seq}\n")
    return out_path

# ---------------- CLI ----------------
def _cli():
    ap = argparse.ArgumentParser(description="Bio files processor")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("convert-fasta", help="FASTA multiline → one-line")
    p1.add_argument("--input-fasta", required=True)
    p1.add_argument("--output-fasta")

    p2 = sub.add_parser("parse-blast", help="BLAST txt → best-hit descriptions (sorted)")
    p2.add_argument("--input-file", required=True)
    p2.add_argument("--output-file", required=True)

    p3 = sub.add_parser("gbk-neighbors", help="Select neighbor CDS translations to FASTA")
    p3.add_argument("--input-gbk", required=True)
    p3.add_argument("--genes", required=True, help="Comma/space-separated string or multiple args", nargs="+")
    p3.add_argument("--n-before", type=int, default=1)
    p3.add_argument("--n-after", type=int, default=1)
    p3.add_argument("--output-fasta", default="neighbors.fasta")

    args = ap.parse_args()
    if args.cmd == "convert-fasta":
        out = convert_multiline_fasta_to_oneline(args.input_fasta, args.output_fasta)
        print(out)
    elif args.cmd == "parse-blast":
        out = parse_blast_output(args.input_file, args.output_file)
        print(out)
    elif args.cmd == "gbk-neighbors":
        # accept either a list of tokens or a single long string
        genes = args.genes if len(args.genes) > 1 else args.genes[0]
        out = select_genes_from_gbk_to_fasta(
            input_gbk=args.input_gbk,
            genes=genes,
            n_before=args.n_before,
            n_after=args.n_after,
            output_fasta=args.output_fasta,
        )
        print(out)

if __name__ == "__main__":
    _cli()
