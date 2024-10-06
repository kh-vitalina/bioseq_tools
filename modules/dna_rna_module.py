def transcribe(sequence: str) -> str:
    return sequence.replace("T", "U").replace("t", "u")


def reverse(sequence: str) -> str:
    return sequence[::-1]


def complement(sequence: str) -> str:
    complement_dict = {
        "A": "T", "T": "A", "C": "G", "G": "C", "U": "A",
        "a": "t", "t": "a", "c": "g", "g": "c", "u": "a",
    }
    return "".join([complement_dict[nuc] for nuc in sequence])


def reverse_complement(sequence: str) -> str:
    return reverse(complement(sequence))
