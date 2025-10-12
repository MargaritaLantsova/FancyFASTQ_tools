<div align="center">
    
# FancyFASTQ_tools

FancyFASTQ_tools is a comprehensive Python toolkit for bioinformatics sequence analysis, providing utilities for DNA/RNA sequence manipulation and FASTQ file filtering. The package offers both programmatic and interactive interfaces for bioinformatics workflows.

![Bioinformatics](https://wallpapers.com/images/hd/dna-background-np9ddwjkk0tcwe36.jpg)
*Источник: [wallpapers.com](https://wallpapers.com/background/dna-background-np9ddwjkk0tcwe36.html)*

</div>
## Content

- [Installation](#installation)
- [Features](#features)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [API Reference](#api-reference)
- [FAQ](#faq)
- [Project Structure](#project-structure)


## Project Structure

```
FancyFASTQ_tools/
├── main.py
├── README.md
└── modules/
    ├── __init__.py
    ├── sequence_tools.py
    └── fastq_tools.py
```

**File descriptions:**
- `main.py` - Main entry point and command-line interface
- `README.md` - Project documentation
- `modules/__init__.py` - Package initialization
- `modules/sequence_tools.py` - DNA/RNA sequence operations
- `modules/fastq_tools.py` - FASTQ processing utilities

## Installation
### Prerequisites
- Python 3.6 or higher
- No external dependencies required (uses only standard library)

### Installation from GitHub

```bash
git clone https://github.com/your_username/FancyFASTQ_tools.git
cd FancyFASTQ_tools
```

## Features
## DNA/RNA Sequence Tools

- 🔄 **Reverse**: Generate reverse sequences
- 🧬 **Complement**: Generate complementary sequences
- ⚡ **Reverse Complement**: Generate reverse-complement sequences
- 📝 **Transcription**: Convert between DNA and RNA
- ✅ **Sequence Validation**: Validate nucleic acid sequences
- 📊 **GC Content**: Calculate GC composition percentage
- 🔍 **Type Detection**: Identify DNA vs RNA sequences

## FASTQ Processing

- 🎯 **Quality Filtering**: Filter by average Phred quality scores
- 📈 **GC Content Filtering**: Filter sequences by GC composition
- 📏 **Length Filtering**: Filter sequences by length boundaries
- 🎛️ **Multi-criteria Filtering**: Combine multiple filtering criteria


## Quick Start
## Running the Main Application

```bash
# Run the interactive application
python main.py
```

## Interactive Mode
```python
from main import run_dna_rna_tools

# Launch interactive menu system
run_dna_rna_tools()
```
## Programmatic Mode
## DNA/RNA Sequence Manipulation
```python
from main import run_dna_rna_tools, filter_fastq

# DNA/RNA operations
result = run_dna_rna_tools("ATCG", "reverse")  # Returns "GCTA"
result = run_dna_rna_tools("ATCG", "gc_content")  # Returns 50.0

# Multiple sequences
results = run_dna_rna_tools("ATCG", "GCTA", "reverse")  # Returns ["GCTA", "ATCG"]
```
## FASTQ Sequence Filtering
```python
from main import run_dna_rna_tools

# Basic operations
sequence = "ATCGATCG"
reverse_seq = run_dna_rna_tools(sequence, "reverse")
complement_seq = run_dna_rna_tools(sequence, "complement") 
reverse_comp_seq = run_dna_rna_tools(sequence, "reverse_complement")
transcribed = run_dna_rna_tools(sequence, "transcribe")

print(f"Original: {sequence}")
print(f"Reverse: {reverse_seq}")
print(f"Complement: {complement_seq}")
print(f"Reverse Complement: {reverse_comp_seq}")
print(f"Transcribed: {transcribed}")
```
## Usage Examples
## DNA/RNA Sequence Manipulation
```python
from main import run_dna_rna_tools

# Basic operations
sequence = "ATCGATCG"
reverse_seq = run_dna_rna_tools(sequence, "reverse")
complement_seq = run_dna_rna_tools(sequence, "complement") 
reverse_comp_seq = run_dna_rna_tools(sequence, "reverse_complement")
transcribed = run_dna_rna_tools(sequence, "transcribe")

print(f"Original: {sequence}")
print(f"Reverse: {reverse_seq}")
print(f"Complement: {complement_seq}")
print(f"Reverse Complement: {reverse_comp_seq}")
print(f"Transcribed: {transcribed}")
```

## FASTQ Sequence Filtering
```python
from main import filter_fastq

# Example FASTQ data
sequences = {
    "read1": ("ATCGATCGATCG", "IIIIIIIIIIII"),  # GC=50%, quality=40
    "read2": ("ATATATAT", "IIIIIIII"),          # GC=0%, quality=40  
    "read3": ("GGGCCC", "AAAAAA"),              # GC=100%, quality=0
    "read4": ("A" * 100, "I" * 100)            # Length=100, GC=0%
}

# Filter by GC content (30-70%) and quality (≥20)
filtered = filter_fastq(
    sequences, 
    gc_bounds=(30, 70),
    quality_threshold=20
)

print(f"Original: {len(sequences)} sequences")
print(f"Filtered: {len(filtered)} sequences")
```

### Streaming FASTQ filtering (CLI)

Streaming filter: read a FASTQ record → check filters → write immediately. Keeps RAM usage low.

```bash
python main.py fastq-filter \
  --input-fastq path/to/reads.fastq.gz \
  --output-fastq filtered_reads.fastq \
  --gc-bounds 40 65 \
  --length-bounds 50 150 \
  --min-qual 20
```

**Arguments**

* `--input-fastq` — path to input FASTQ (`.fastq`/`.fq` and `.fastq.gz`/`.fq.gz` supported).
* `--output-fastq` — output FASTQ **written to `./filtered/`**.
* `--gc-bounds` — GC bounds: one upper threshold (`60`) or two values (`min max`, e.g. `40 65`).
* `--length-bounds` — length bounds: one upper threshold (`150`) or two values (`50 150`).
* `--min-qual` — minimal average Phred+33.

**Example output**

```
Input: path/to/reads.fastq.gz
Output: filtered/filtered_reads.fastq
Total: 123456
Kept: 78910
```

---

### `bio_files_processor.py` (CLI)

Three subcommands.

#### `convert-fasta` — multiline FASTA → one-line FASTA

```bash
python modules/bio_files_processor.py convert-fasta \
  --input-fasta data/genes.fasta \
  --output-fasta genes.oneline.fasta
```

If `--output-fasta` is omitted, the file is created next to the input with the `.oneline.fasta` suffix.

#### `parse-blast` — extract best-hit descriptions from BLAST TXT

For every section `Sequences producing significant alignments:` take the **first** table row (best hit) and keep the **Description** (first column). Output is a unique, alphabetically sorted list (one per line).

```bash
python modules/bio_files_processor.py parse-blast \
  --input-file results/example_blast_results.txt \
  --output-file best_hits.txt
```

#### `gbk-neighbors` — neighbor CDS translations from GBK → FASTA

From a GBK annotation, for each gene of interest (match by `/gene` or `/locus_tag`) write **n_before** genes before and **n_after** after (exclude the targets themselves).

```bash
python modules/bio_files_processor.py gbk-neighbors \
  --input-gbk data/ecoli.gbk \
  --genes acrA tolC marA \
  --n-before 2 \
  --n-after 2 \
  --output-fasta neighbors.fasta
```

You may also pass a single string with separators:

```bash
python modules/bio_files_processor.py gbk-neighbors \
  --input-gbk data/ecoli.gbk \
  --genes "acrA, tolC; marA" \
  --output-fasta neighbors.fasta
```

---

### Safe output & atomic writes

* **No overwrites:** if a target path exists, a suffix `__1`, `__2`, … is appended until a free name is found.
* **Atomic writes:** data is written to a temp file and then swapped into place via `os.replace`, preventing half-written files.
* **FASTQ filter output:** always saved under `./filtered/` (created if missing) using the name from `--output-fastq`, with the same uniqueness guarantees.

## API Reference
```
run_dna_rna_tools(*args)
```
Process DNA/RNA sequences with specified operations.

Parameters:
```
    *args: Sequences followed by operation name

    Operations: "reverse", "complement", "reverse_complement", "transcribe", "is_dna", "is_rna", "gc_content", "is_nucleic_acid"
```
Returns: Single result or list of results
filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0)

Filter FASTQ sequences by multiple criteria.

Parameters:
```
    seqs: Dictionary {read_name: (sequence, quality)}

    gc_bounds: GC content bounds (number or tuple)

    length_bounds: Length bounds (number or tuple)

    quality_threshold: Minimum average Phred quality
```
Returns: Filtered dictionary of sequences


## FAQ

### ❓ Can I process RNA sequences with this toolkit?
**✅ Yes!** The toolkit automatically detects and handles both DNA and RNA sequences, including proper transcription between them.

### ❓ What quality encoding does the FASTQ filter use?
**✅** The toolkit uses **Phred+33 encoding**, which is standard for most modern sequencing platforms.

### ❓ Can I filter by multiple criteria simultaneously?
**✅ Absolutely!** The `filter_fastq` function applies all specified filters (GC content, length, and quality) together.

### ❓ Does the toolkit handle mixed-case sequences?
**✅ Yes,** all sequence operations are case-insensitive and preserve the original case in output.

### ❓ What happens if I have both T and U in a sequence?
**⚠️** The validation functions will identify this as an **invalid nucleic acid sequence**.

### ❓ How do I run the application?
**🚀** Simply run `python main.py` from the command line in the FancyFASTQ_tools directory.
