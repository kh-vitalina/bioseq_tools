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

