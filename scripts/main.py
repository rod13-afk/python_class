"""
NAME
       Arguments AT.
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This program uses a text file that is passed as an argument from the command line,
        this file contains DNA sequences, the program opens it and calculates the AT percentages
        of the sequences that do not contain N's. The results are saved in a new text file and
        the sequences that had N's and the amount of these are printed on the terminal.
CATEGORY
        Genomic sequence.
USAGE
        python3 arguments_AT.py -i data/4_dna_sequences.txt -o output/dna_sequences.txt
ARGUMENTS
        -h, --help            show this help message and exit
        -i, --input           the input arguments
        -o, --output          output results.
        --round               use integer values
INPUT
        A fasta file with DNA sequences.
OUTPUT
        It returns as output a text file with the AT percentages of each sequence and prints
        the name of the sequences containing N's and the quantity of these in the terminal.
EXAMPLES
        Input
            python3 arguments_AT.py -i data/4_dna_sequences.txt -o output/dna_sequences.txt
         Output
             A dna_sequences.txt file in the output folder.
"""

import sequenceAnalysis
import argparse

# Create the parser
parser = argparse.ArgumentParser(description="This program calculates the content of AT, skipping "
                                             "lines that contain N's and taking arguments from the command line")

# Add the arguments
parser.add_argument('-i', '--input',
                    metavar="path/to/file",
                    help="File with gene sequences",
                    required=True)
parser.add_argument('-o', '--output',
                    metavar="path/to/output",
                    help="Path for the output file",
                    required=False)

# Execute the parse_args method
arguments = parser.parse_args()

# Code that calls the functions
file = sequenceAnalysis.open_file(arguments.input)
rna = sequenceAnalysis.transcription(file)
sequenceAnalysis.get_AT_content(file)
sequenceAnalysis.create_fasta(arguments.output, rna)
