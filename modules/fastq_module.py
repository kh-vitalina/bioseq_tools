def count_gc(sequence: str) -> float:
    g_count = sequence.count("G") + sequence.count("g")
    c_count = sequence.count("C") + sequence.count("c")
    return (g_count + c_count) / len(sequence) * 100


def get_mean_quality(quality: str) -> float:
    return sum(ord(q) - 33 for q in quality) / len(quality)


def read_fastq(file):
    while True:
        seq_id = file.readline().strip()
        seq = file.readline().strip()
        file.readline()  #skip line with +
        quality = file.readline().strip()
        if not seq_id:
            break
        yield seq_id, seq, quality


def write_fastq(file, seq_id, seq, quality):
    file.write(f"{seq_id}\n")
    file.write(f"{seq}\n")
    file.write("+\n")
    file.write(f"{quality}\n")
