"""
This package contains core modules for bioinformatics sequence analysis:
- sequence_tools: DNA/RNA sequence manipulation utilities
- fastq_tools: FASTQ file processing and filtering utilities
"""

from .sequence_tools import (
    is_valid_nucleic_acid,
    is_dna,
    is_rna,
    transcribe,
    reverse,
    complement,
    reverse_complement,
    gc_content,
    OPERATIONS,
)

from .fastq_tools import (
    calculate_avg_quality,
    check_bounds,
    filter_sequence,
    DEFAULT_GC_BOUNDS,
    DEFAULT_LENGTH_BOUNDS,
    DEFAULT_QUALITY_THRESHOLD,
)

# Define what gets imported with "from modules import *"
__all__ = [
    # Sequence tools
    "is_valid_nucleic_acid",
    "is_dna",
    "is_rna",
    "transcribe",
    "reverse",
    "complement",
    "reverse_complement",
    "gc_content",
    "OPERATIONS",
    # FASTQ tools
    "calculate_avg_quality",
    "check_bounds",
    "filter_sequence",
    "DEFAULT_GC_BOUNDS",
    "DEFAULT_LENGTH_BOUNDS",
    "DEFAULT_QUALITY_THRESHOLD",
]
