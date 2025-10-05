def is_valid_nucleic_acid(sequence):
    """Check if sequence is valid DNA or RNA"""
    if not sequence:
        return False

    upper_seq = sequence.upper()
    valid_dna = set("ATCG")
    valid_rna = set("AUCG")
    seq_set = set(upper_seq)

    # Check if sequence contains only valid nucleotides
    if not (seq_set.issubset(valid_dna) or seq_set.issubset(valid_rna)):
        return False

    # Check for mixed T/U
    has_t = "T" in upper_seq
    has_u = "U" in upper_seq

    return not (has_t and has_u)


def is_dna(sequence):
    """Check if sequence is DNA"""
    if not is_valid_nucleic_acid(sequence):
        return False

    upper_seq = sequence.upper()
    return "T" in upper_seq or (
        "U" not in upper_seq and set(upper_seq).issubset(set("ATCG"))
    )


def is_rna(sequence):
    """Check if sequence is RNA"""
    if not is_valid_nucleic_acid(sequence):
        return False

    upper_seq = sequence.upper()
    return "U" in upper_seq or (
        "T" not in upper_seq and set(upper_seq).issubset(set("AUCG"))
    )


def transcribe(sequence):
    """Transcribe DNA↔RNA"""
    if not is_valid_nucleic_acid(sequence):
        raise ValueError(f"Invalid nucleic acid sequence: {sequence}")

    transcription_map = {
        "T": "U",
        "t": "u",
        "U": "T",
        "u": "t",
        "A": "A",
        "a": "a",
        "C": "C",
        "c": "c",
        "G": "G",
        "g": "g",
    }

    return "".join(transcription_map[nucleotide] for nucleotide in sequence)


def reverse(sequence):
    """Return reversed sequence"""
    return sequence[::-1]


def complement(sequence):
    """Return complementary sequence"""
    if not is_valid_nucleic_acid(sequence):
        raise ValueError(f"Invalid nucleic acid sequence: {sequence}")

    is_dna_seq = is_dna(sequence)

    dna_complement_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }

    rna_complement_map = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }

    complement_map = dna_complement_map if is_dna_seq else rna_complement_map
    return "".join(complement_map.get(nuc, nuc) for nuc in sequence)


def reverse_complement(sequence):
    """Return reverse complement sequence"""
    return complement(reverse(sequence))


def gc_content(sequence):
    """Calculate GC content in percentage"""
    if not sequence:
        return 0.0

    upper_seq = sequence.upper()
    gc_count = upper_seq.count("G") + upper_seq.count("C")
    return (gc_count / len(upper_seq)) * 100


# Dictionary of available operations
OPERATIONS = {
    "is_nucleic_acid": lambda seq: True,
    "transcribe": transcribe,
    "reverse": reverse,
    "complement": complement,
    "reverse_complement": reverse_complement,
    "is_dna": is_dna,
    "is_rna": is_rna,
    "gc_content": gc_content,
}
