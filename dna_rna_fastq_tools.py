def run_dna_rna_tools(*args):
    procedure = args[-1]
    sequences = args[:-1]

    for seq in sequences:
        if not set(seq).issubset(set("ATCGUatcgu")) or ('T' in seq and 'U' in seq):
            return f"Uncorrect sequence (not DNA or RNA): {seq}"

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


def transcribe(sequence):
    return sequence.replace("T", "U").replace("t", "u")


def reverse(sequence):
    return sequence[::-1]


def complement(sequence):
    complement_dict = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
        "a": "t",
        "t": "a",
        "c": "g",
        "g": "c",
        "u": "a",
    }
    return "".join([complement_dict[nuc] for nuc in sequence])


def reverse_complement(sequence):
    return reverse(complement(sequence))


def count_gc(sequence):
    g_count = sequence.count('G') + sequence.count('g')
    c_count = sequence.count('C') + sequence.count('c')
    return (g_count + c_count) / len(sequence) * 100


def get_mean_quality(quality):
    return sum(ord(q) - 33 for q in quality) / len(quality)


def filter_fastq(seqs, gc_bounds = (0, 100), length_bounds = (0, 2 ** 32), quality_threshold = 0):
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
        
    gc_min, gc_max = gc_bounds
    length_min, length_max = length_bounds

    filtered_seqs = {}
    
    for seq_id, (seq, quality) in seqs.items():
        gc_content = count_gc(seq)
        if gc_min <= gc_content <= gc_max and length_min <= len(seq) <= length_max and get_mean_quality(quality) >= quality_threshold:
            filtered_seqs[seq_id] = (seq, quality)
            
    return filtered_seqs
