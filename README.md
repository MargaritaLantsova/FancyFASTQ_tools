<div align="center">

# FancyFASTQ_tools

A Python toolkit for bioinformatics sequence analysis built around:

* 🧬 Object-Oriented DNA/RNA/Protein sequence classes
* 📈 FASTQ filtering using **Biopython**
* 📂 Utility tools for FASTA, BLAST and GenBank processing

</div>

---

# Content

* [Installation](#installation)
* [Project Structure](#project-structure)
* [Sequence Classes (OOP)](#sequence-classes-oop)
* [FASTQ Filtering (Biopython)](#fastq-filtering-biopython)
* [File Processing Utilities](#file-processing-utilities)
* [Command Line Usage](#command-line-usage)

---

# Project Structure

```
FancyFASTQ_tools/
├── main.py
├── README.md
└── requirements.txt
```

**File descriptions:**

* `main.py` — all sequence classes, FASTQ filtering and CLI
* `requirements.txt` — external dependencies

There are no legacy modules. All functionality is consolidated into `main.py`.



# Installation

## Requirements

* Python 3.8+
* Biopython

## Install dependencies

```bash
pip install -r requirements.txt
```

`requirements.txt` contains:

```
biopython>=1.80
```

# Sequence Classes (OOP)

The toolkit now uses proper object-oriented design.

## Available classes

* `BiologicalSequence` (abstract)
* `NucleicAcidSequence` (abstract)
* `DNASequence`
* `RNASequence`
* `AminoAcidSequence`

## Example usage

```python
from main import DNASequence, RNASequence, AminoAcidSequence

dna = DNASequence("ATGCGT")
print(dna.reverse())
print(dna.complement())
print(dna.reverse_complement())
print(dna.gc_content())

rna = dna.transcribe()
print(rna)

protein = AminoAcidSequence("MKWVTFISLL")
print(protein.molecular_weight_approx())
```

### Features

DNA / RNA:

* Reverse
* Complement
* Reverse complement
* Transcription (DNA → RNA)
* GC content
* Alphabet validation
* Slicing returns same object type

Protein:

* Alphabet validation
* Approximate molecular weight calculation

# FASTQ Filtering (Biopython)

FASTQ filtering is implemented using **Biopython (SeqIO, SeqUtils)**.

Filters applied:

* Length bounds
* Mean Phred quality
* GC content percentage

## CLI Example

```bash
python main.py fastq-filter \
  --input-fastq reads.fastq.gz \
  --output-fastq filtered.fastq \
  --gc-bounds 40 65 \
  --length-bounds 50 150 \
  --min-qual 20
```

### Arguments

* `--input-fastq` — input `.fastq` or `.fastq.gz`
* `--output-fastq` — output filename (saved to `./filtered/`)
* `--gc-bounds` — one value (upper) or two values (`min max`)
* `--length-bounds` — one value or two values
* `--min-qual` — minimum mean Phred score

### Example Output

```
Input: reads.fastq.gz
Output: filtered/filtered.fastq
Total: 120000
Kept: 84213
```

# File Processing Utilities

All utilities are integrated into the main CLI.

## 1. Convert FASTA (multiline → one-line)

```bash
python main.py convert-fasta \
  --input-fasta genes.fasta \
  --output-fasta genes.oneline.fasta
```

## 2. Parse BLAST output (best hits)

Extract the first hit from each:

```
Sequences producing significant alignments:
```

section and output unique descriptions (sorted).

```bash
python main.py parse-blast \
  --input-file blast_results.txt \
  --output-file best_hits.txt
```


## 3. GBK neighbor CDS extraction

For each gene of interest, extract neighboring CDS translations.

```bash
python main.py gbk-neighbors \
  --input-gbk ecoli.gbk \
  --genes acrA tolC marA \
  --n-before 2 \
  --n-after 2 \
  --output-fasta neighbors.fasta
```

You may also pass:

```bash
--genes "acrA, tolC; marA"
```

# Command Line Overview

Run help:

```bash
python main.py --help
```

Available commands:

* `fastq-filter`
* `convert-fasta`
* `parse-blast`
* `gbk-neighbors`


# Notes

* No legacy modules remain.
* FASTQ filtering relies on Biopython.
* All sequence logic is implemented via OOP with polymorphism.
* Slicing preserves sequence type (e.g. DNA → DNA).
* Validation raises exceptions on invalid symbols.

---
