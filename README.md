# Bioinformatics Tools

This package provides utilities for DNA/RNA sequence processing and filtering of FASTQ sequences, and processing of FASTA and BLAST output data.

## Installation

Clone the repository and navigate to the project directory:

`bash`
```sh
git clone https://github.com/kh-vitalina/bioseq_tools.git && cd bioseq_tools
```

## Usage:

## dna_rna_fastq_tools.py

**DNA/RNA Processing** 
- The main function run_dna_rna_tools allows you to perform various DNA/RNA operations.

#### Example:

~~~python
from dna_rna_fastq_tools import run_dna_rna_tools
result = run_dna_rna_tools("AtGcTCgGtA", "complement")
print(result) # Output: "TaCgAGcCaT"
~~~

**FASTQ Filtering**
- The function filter_fastq filters sequences from a FASTQ file based on specified criteria.

#### Example:

~~~python
from dna_rna_fastq_tools import filter_fastq

input_fastq = "input_sequences.fastq"
output_fastq = "filtered_sequences.fastq"

filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), quality_threshold=20)
~~~

#### Overview:

The package includes two main functionalities:
- **run_dna_rna_tools**: Functions for transcribing, reversing, complementing, and getting reverse complements of sequences.
- **filter_fastq**: Filters sequences based on GC bounds, length, and quality threshold.

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

## Modules Structure

#### modules
- **dna_rna_module.py**: contains functions related to DNA/RNA processing.
- **fastq_module.py**: contains functions for filtering FASTQ sequences.
- **read_fastq**: generator function to read sequences from a FASTQ file
- **write_fastq**: write sequences to a FASTQ file
