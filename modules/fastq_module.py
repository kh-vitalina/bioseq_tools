def count_gc(sequence: str) -> float:
    g_count = sequence.count("G") + sequence.count("g")
    c_count = sequence.count("C") + sequence.count("c")
    return (g_count + c_count) / len(sequence) * 100


def get_mean_quality(quality: str) -> float:
    return sum(ord(q) - 33 for q in quality) / len(quality)
