from .sequence_tools import gc_content

# Определяем константы прямо в модуле
DEFAULT_GC_BOUNDS = (0, 100)
DEFAULT_LENGTH_BOUNDS = (0, 2**32)
DEFAULT_QUALITY_THRESHOLD = 0


def calculate_avg_quality(quality_string: str) -> float:
    """
    Calculate average read quality in phred33
    quality_string: string with quality characters
    """
    if not quality_string:
        return 0.0
    
    quality_scores = [ord(char) - 33 for char in quality_string]
    return sum(quality_scores) / len(quality_scores)


def check_bounds(value: float, bounds) -> bool:
    """
    Check if value is within specified bounds
    bounds can be tuple (min, max) or single number (max)
    """
    if isinstance(bounds, (int, float)):
        return 0 <= value <= bounds
    elif isinstance(bounds, tuple) and len(bounds) == 2:
        return bounds[0] <= value <= bounds[1]
    else:
        raise ValueError("Bounds must be a number or tuple of two numbers")


def filter_sequence(seq_data, gc_bounds=DEFAULT_GC_BOUNDS, 
                   length_bounds=DEFAULT_LENGTH_BOUNDS, 
                   quality_threshold=DEFAULT_QUALITY_THRESHOLD):
    """
    Filter single sequence by specified criteria
    seq_data: tuple (sequence, quality)
    """
    sequence, quality = seq_data
    
    # Check length
    seq_length = len(sequence)
    if not check_bounds(seq_length, length_bounds):
        return False
    
    # Check GC content
    gc_percent = gc_content(sequence)
    if not check_bounds(gc_percent, gc_bounds):
        return False
    
    # Check quality
    avg_quality = calculate_avg_quality(quality)
    if avg_quality < quality_threshold:
        return False
    
    return True


