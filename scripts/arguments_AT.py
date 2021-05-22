'''
NAME
       Arguments AT.
VERSION
        [1.0]
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
        A text file with DNA sequences.
OUTPUT
        It returns as output a text file with the AT percentages of each sequence and prints
        the name of the sequences containing N's and the quantity of these in the terminal.
EXAMPLES
        Input
            python3 arguments_AT.py -i data/4_dna_sequences.txt -o output/dna_sequences.txt
         Output
             A dna_sequences.txt file in the output folder.
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/scripts/arguments_AT.py

'''

# The library needed to work with arguments is imported
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

# Print arguments
print(arguments.input, arguments.output)


# An AmbiguousBaseError class is created to avoid hiding errors.
class AmbiguousBaseError(Exception):
    pass


# Function to calculate the at content
def get_at_content(sequence, default=2):
    if sequence.count("N") > 0:
        raise AmbiguousBaseError("Sequence contains" + str(sequence.count("N")) + "N's")

    length = len(sequence) - 1
    a_content = (sequence.upper()).count("A")
    t_content = (sequence.upper()).count("T")
    at_content = (a_content + t_content) / length
    return round(at_content, default)


# A 'try' to open the requested file. If the file is not there, the user is asked to enter the name of the file
try:
    input_file = open(arguments.input)
except IOError as ex:
    print("Sorry, couldn't find the file: " + ex.strerror)
    alternative_file = input("Please enter your file: ")
    input_file = open(alternative_file)

file_to_work = input_file.readlines()
input_file.close()


# The output file is created in the output folder.
output_file = open("output/4_AT_content.txt", "w")

for i, line in enumerate(file_to_work):

    # Each line is divided.
    line = line.split(' = ')

    # The hyphens and quotes are removed from the DNA sequences
    non_dash = str(line[1]).replace('-', '')
    non_quotes = non_dash.replace('"', '')

    # A try to save the new message of each sequence in the output. If the sequence contains N's,
    # it is skipped and the number of N's is printed.
    try:
        new_message = "SEQ_" + str(i + 1) + " has " + str(get_at_content(non_quotes)) + " AT content" + "\n"
        output_file.write(new_message)
    except AmbiguousBaseError as ex:
        print("Skipping invalid sequence SEQ_" + str(i + 1) + ex.args[0])

# The output file is closed
output_file.close()
