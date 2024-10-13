from modules.dna_rna_module import transcribe, reverse, complement, reverse_complement
from modules.fastq_module import count_gc, get_mean_quality, read_fastq, write_fastq


def run_dna_rna_tools(*args: str) -> str:
    """Runs DNA/RNA processing procedures on prescribed sequences"""
    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not set(seq).issubset(set("ATCGUatcgu")) or ("T" in seq and "U" in seq):
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


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: float = 0,
):
    """
    Filters sequences from input FASTQ file based on GC bounds, length and quality threshold
    Args:
        input_fastq: path to the input FASTQ file
        output_fastq: path to save filtered sequences
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, int):
        length_bounds = (length_bounds, length_bounds)

    gc_min, gc_max = gc_bounds
    length_min, length_max = length_bounds

    with open(input_fastq, "r") as infile, open(
        f"filtered/{output_fastq}", "w"
    ) as outfile:
        for record in read_fastq(infile):
            seq_id, seq, quality = record
            if (
                gc_min <= count_gc(seq) <= gc_max
                and length_min <= len(seq) <= length_max
                and get_mean_quality(quality) >= quality_threshold
            ):
                write_fastq(outfile, seq_id, seq, quality)
