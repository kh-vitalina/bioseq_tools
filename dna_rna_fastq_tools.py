from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os

class BiologicalSequence(ABC):
    """Base class for all biological sequences"""
    def __init__(self, sequence):
        self.sequence = sequence.upper()
        self._check_alphabet()
    
    @abstractmethod
    def _check_alphabet(self):
        """Verify that sequence doesn't contain invalid characters"""
        pass
    
    def __len__(self):
        return len(self.sequence)
    
    def __getitem__(self, index):
        """Get one character or a slice of the sequence"""
        result = self.sequence[index]
        if isinstance(index, slice):
            return self.__class__(result)
        return result
    
    def __str__(self):
        return self.sequence
    
    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"


class NucleicAcidSequence(BiologicalSequence):
    """Base class for DNA and RNA sequences"""
    complement_table = {}
    
    def _check_alphabet(self):
        """Validate nucleic acid alphabet"""
        valid = set(self.complement_table.keys())
        for letter in self.sequence:
            if letter not in valid:
                raise ValueError(f"Invalid letter '{letter}' in nucleic acid sequence")
    
    def complement(self):
        """Return complementary sequence"""
        if not self.complement_table:
            raise NotImplementedError("use DNASequence or RNASequence")
        
        complemented = ''.join(self.complement_table[base] for base in self.sequence)
        return self.__class__(complemented)
    
    def reverse(self):
        """Return reversed sequence"""
        return self.__class__(self.sequence[::-1])
    
    def reverse_complement(self):
        """Return reverse complement"""
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """DNA sequence class"""  
    complement_table = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'
    }
    
    def transcribe(self):
        """Transcribe DNA to RNA"""
        rna = self.sequence.replace('T', 'U')
        return RNASequence(rna)


class RNASequence(NucleicAcidSequence):
    """RNA sequence class"""  
    complement_table = {
        'A': 'U', 'U': 'A',
        'G': 'C', 'C': 'G',
        'N': 'N'
    }


class AminoAcidSequence(BiologicalSequence):
    """Protein sequence class"""
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    def _check_alphabet(self):
        """Validate amino acid alphabet"""
        for letter in self.sequence:
            if letter not in self.valid_amino_acids:
                raise ValueError(f"Invalid amino acid '{letter}' in protein sequence")
    
    def molecular_weight(self):
        """Calculate molecular weight"""
        #molecular weights g/mol
        weights = {
            'A': 89, 'C': 121, 'D': 133, 'E': 147,
            'F': 165, 'G': 75, 'H': 155, 'I': 131,
            'K': 146, 'L': 131, 'M': 149, 'N': 132,
            'P': 115, 'Q': 146, 'R': 174, 'S': 105,
            'T': 119, 'V': 117, 'W': 204, 'Y': 181
        }
        total = sum(weights.get(aa, 0) for aa in self.sequence)
        total -= (len(self.sequence) - 1) * 18
        return total
    

def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), 
                 length_bounds=(0, 2**32), quality_threshold=0):
    """
    Filter FASTQ reads by GC% content, length and quality
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (length_bounds, length_bounds)
    gc_min, gc_max = gc_bounds
    len_min, len_max = length_bounds
    
    output_dir = os.path.dirname(output_fastq)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    filtered = 0
    total = 0
    with open(output_fastq, 'w') as out:
        for record in SeqIO.parse(input_fastq, 'fastq'):
            total += 1
            seq = str(record.seq)
            if len(seq) < len_min or len(seq) > len_max:
                continue
            gc_content = gc_fraction(record.seq) * 100
            if gc_content < gc_min or gc_content > gc_max:
                continue
            if quality_threshold > 0:
                qualities = record.letter_annotations.get('phred_quality', [])
                if qualities:
                    avg_quality = sum(qualities) / len(qualities)
                    if avg_quality < quality_threshold:
                        continue
            
            SeqIO.write(record, out, 'fastq')
            filtered += 1
    
    print(f"Processed {total} reads, kept {filtered}")
    return filtered