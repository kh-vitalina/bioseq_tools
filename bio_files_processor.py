def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):

    """
    Converts a multiline FASTA file to oneline FASTA file
    
    Args:
        input_fasta: path to the input FASTA file
        output_fasta: path to the output FASTA file (optional. If not provide, output file will have the same name as the input file with addition “_oneline.fasta”)

    Return:
        None
    """
    
    if output_fasta is None:
        output_fasta = f"{input_fasta.rsplit('.', 1)[0]}_oneline.fasta"
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence = ''
        header = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    outfile.write(f"{header}\n{sequence}\n")
                header = line
                sequence = ''
            else:
                sequence += line
        if sequence:
            outfile.write(f"{header}\n{sequence}\n")


def parse_blast_output(input_file: str, output_file: str):
    """
    Parse BLAST output file and extract the best match description for each query

    Args:
        input_file: path to the BLAST output file (txt format)
        output_file: path to the output file where to save results 

    Returns:
        None
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        best_match = []
        in_significant_alignments = False
        
        for line in infile:
            if "Sequences producing significant alignments:" in line:
                in_significant_alignments = True
                continue
            if in_significant_alignments:
                if line.strip():
                    columns = line.split()
                    if columns:
                        best_match.append(columns[0].strip())

        proteins = sorted(set(best_match))
        
        for protein in proteins:
            outfile.write(f"{protein}\n")
