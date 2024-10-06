from modules.dna_rna_module import transcribe, reverse, complement, reverse_complement
from modules.fastq_module import count_gc, get_mean_quality


def run_dna_rna_tools(*args: str) -> str:
    """Runs DNA/RNA processing procedures on prescribed sequences"""
    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not set(seq).issubset(set("ATCGUatcgu")) or ('T' in seq and 'U' in seq):
            return f"Incorrect sequence (not DNA or RNA): {seq}"

    if procedure == "transcribe":
        fn = transcribe
    elif procedure == "reverse":
        fn = reverse
    elif procedure == "complement":
        fn = complement
    elif procedure == "reverse_complement":
        fn = reverse_complement
    else:
        return f"Unknown procedure: {procedure}"

    result = [fn(seq) for seq in sequences]

    if len(sequences) == 1:
        return result[0]
    return result


def filter_fastq(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2 ** 32), quality_threshold: float = 0) -> dict:
    """Filters FASTQ sequences based on GC bounds, length and quality threshold"""
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    gc_min, gc_max = gc_bounds
    length_min, length_max = length_bounds

    filtered_seqs = {}

    for seq_id, (seq, quality) in seqs.items():
        gc_content = count_gc(seq)
        if (gc_min <= gc_content <= gc_max and length_min <= len(seq) <= length_max and get_mean_quality(quality) >= quality_threshold):
            filtered_seqs[seq_id] = (seq, quality)

    return filtered_seqs
