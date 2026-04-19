import pytest
import os
import tempfile
from dna_rna_fastq_tools import (
    DNASequence, RNASequence, AminoAcidSequence, filter_fastq
)


def make_fastq_file(filepath, records):
    """Helper to create a test FASTQ file"""
    with open(filepath, 'w') as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")


class TestDNASequence:

    def test_complement(self):
        """Check that complement works correctly"""
        dna = DNASequence("ATGC")
        assert str(dna.complement()) == "TACG"

    def test_transcribe(self):
        """Check DNA to RNA transcription"""
        dna = DNASequence("ATGC")
        rna = dna.transcribe()
        assert str(rna) == "AUGC"

    def test_invalid_letter(self):
        """Wrong letter raises an error"""
        with pytest.raises(ValueError):
            DNASequence("ATGCX")


class TestAminoAcidSequence:

    def test_molecular_weight(self):
        """Weight of one amino acid"""
        protein = AminoAcidSequence("A")
        assert protein.molecular_weight() == 89

    def test_invalid_amino_acid(self):
        """Invalid amino acid raises an error"""
        with pytest.raises(ValueError):
            AminoAcidSequence("MVXLS")


class TestFilterFastq:

    def setup_method(self):
        self.test_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.test_dir, "input.fastq")
        self.output_file = os.path.join(self.test_dir, "output.fastq")

        #read1 all GC, high quality (I= phred 40)
        #read2 no GC at all, high quality
        #read3 low GC, bad quality (!= phred 0)
        records = [
            ("read1", "GCGCGC", "IIIIII"),
            ("read2", "ATATAT", "IIIIII"),
            ("read3", "GCATAT", "!!!!!!"),
        ]
        make_fastq_file(self.input_file, records)

    def teardown_method(self):
        for f in [self.input_file, self.output_file]:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(self.test_dir):
            os.rmdir(self.test_dir)

    def test_no_filter(self):
        """Without any filters all reads will pass"""
        result = filter_fastq(self.input_file, self.output_file)
        assert result == 3

    def test_gc_filter(self):
        """Only reads with 30-50% GC content stay"""
        result = filter_fastq(self.input_file, self.output_file,
                              gc_bounds=(30, 50))
        assert result == 1

    def test_quality_filter(self):
        """Filtered out low quality read"""
        result = filter_fastq(self.input_file, self.output_file,
                              quality_threshold=20)
        assert result == 2

    def test_length_filter(self):
        """All reads are length 6, so with min 10 none pass"""
        result = filter_fastq(self.input_file, self.output_file,
                              length_bounds=(10, 100))
        assert result == 0

    def test_output_file_exists(self):
        """Check that output file is created and non-empty"""
        filter_fastq(self.input_file, self.output_file)
        assert os.path.exists(self.output_file)
        assert os.path.getsize(self.output_file) > 0

    def test_file_not_found(self):
        """Error if input file does not exist"""
        with pytest.raises(FileNotFoundError):
            filter_fastq("no_such_file.fastq", self.output_file)