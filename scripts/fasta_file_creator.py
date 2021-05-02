'''
NAME
       Programa que crea un archivo fasta.
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This program uses the directory of a text file whose content is DNA sequences,
        opens it and outputs a file with the sequences in the corresponding format.
CATEGORY
        Genomic sequence.
INPUT
        A text file with DNA sequences.
OUTPUT
        Returns as output a fasta file with the sequences in the corresponding format.
EXAMPLES
        Input
            dna/4_dna_sequences.txt
         Output
             A dna_sequences.fasta file in the output folder.
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/scripts/fasta_file_creator.py

'''


# Get the file directory
input_file_directory = input("Insert the directory where your file is: ")

# Here the extension of the file is validated, it has to be .txt if not, the program stops.
validation = input_file_directory.find(".txt")
if validation == -1:
    print("The input file is invalid.")
    exit(1)

# If the input is .txt, the file opens
input_file = open(input_file_directory)

# Here all the lines of the file are stored in the variable all_lines
all_lines = input_file.readlines()

# The file is closed
input_file.close()

# The fasta file is created in the output folder
fasta_file = open("output/dna_sequences.fasta", "x")

''' For each line a head and the sequence are added in the output file.
In this for 'enumerate' is implemented which gives us the opportunity
to count the for loops and use this value in other functionalities of the program'''
for i, line in enumerate(all_lines):

    # Each line is splitted, the first one is the head of the DNA sequence
    line = line.split(' = ')

    # The head of each line is added to the fasta file with a '>' symbol, along with a new line.
    head = "> " + line[0]
    fasta_file.write(head + "\n")

    # This additional functionality tells us if a given sequence has gaps and how many
    if str(line[1]).count('-'):
        print("There are " + str(str(line[1]).count('-')) + " gaps in the " + str(i + 1) + " sequence")

    # Now the hyphens and quotes are removed from the DNA sequences
    non_dash = str(line[1]).replace('-', '')
    non_quotes = non_dash.replace('"', '')

    # The DNA sequences are added to the fasta file.
    fasta_file.write(non_quotes.upper())

# The fasta file is closed
fasta_file.close()

# A message confirming the creation of the file
print("The fasta file was created successfully! Check your output folder.")
