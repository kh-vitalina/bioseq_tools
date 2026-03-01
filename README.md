# Bioinformatics Tools

This package provides utilities for DNA/RNA sequence processing and filtering of FASTQ sequences, and processing of FASTA and BLAST output data.

## Installation

Clone the repository and navigate to the project directory:

`bash`
```sh
git clone https://github.com/kh-vitalina/bioseq_tools.git && cd bioseq_tools
```

## Requirements

```
biopython>=1.79
```

## Usage:

## dna_rna_fastq_tools.py

**DNA/RNA Processing** 
The module now uses classes for working with biological sequences:

- `DNASequence` - for DNA sequences (complement, reverse, reverse_complement, transcribe)
- `RNASequence` - for RNA sequences (complement, reverse, reverse_complement)
- `AminoAcidSequence` - for protein sequences (molecular_weight calculation)

#### Example:

~~~python
from dna_rna_fastq_tools import DNASequence, RNASequence, AminoAcidSequence

#DNA sequence
dna = DNASequence("ATGC")
print(dna.complement())        #output TACG
print(dna.reverse())           #output CGTA
print(dna.reverse_complement()) #output GCAT
print(dna.transcribe())        #output AUGC

#RNA sequence
rna = RNASequence("AUGC")
print(rna.complement())        #output UACG

#Protein sequence
protein = AminoAcidSequence("MVLSPADKT")
print(protein.molecular_weight()) #output molecular weight
~~~

**FASTQ Filtering**
- The function `filter_fastq` filters sequences from a FASTQ file based on specified criteria using Biopython.

#### Example:

~~~python
from dna_rna_fastq_tools import filter_fastq

input_fastq = "input_sequences.fastq"
output_fastq = "filtered_sequences.fastq"

filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), quality_threshold=20)
~~~

#### Overview:

The package includes two main functionalities:
- **Sequence classes**: Object-oriented approach for DNA, RNA and protein sequences
- **filter_fastq**: Filters FASTQ reads based on GC content, length, and quality using Biopython

## bio_files_processor.py

**FASTA Processing**
- The function convert_multiline_fasta_to_oneline converts a multiline FASTA file into a oneline format.

#### Example:

~~~python
from bio_files_processor import convert_multiline_fasta_to_oneline
convert_multiline_fasta_to_oneline("input.fasta", "output_oneline.fasta")
~~~

**BLAST Output Parsing**
- The function parse_blast_output extracts the best match description from a BLAST output file.

#### Example:

~~~python
from bio_files_processor import parse_blast_output
parse_blast_output("blast_output.txt", "parsed_proteins.txt")
~~~

#### Overview:

The package includes two main functionalities:
- **convert_multiline_fasta_to_oneline**: Converts FASTA files with multiline sequences to a oneline format.
- **parse_blast_output**: Parses BLAST output files and extracts best match descriptions.

## Project Structure

- **dna_rna_fastq_tools.py**: main script with sequence classes and FASTQ filtering
- **bio_files_processor.py**: functions for FASTA and BLAST processing
```