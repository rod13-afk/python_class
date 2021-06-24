"""
NAME
       fastaFiles.py
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This module works with fasta files, it contains a function that opens and reads
        fasta files, and its output is a list with the lines of the file. The second
        function creates fasta files, receives as parameters the path and name of the file
        to be created and an iterable object, such as a list, containing the sequences,
        which are the content of the new fasta file.
CATEGORY
        fasta files.
USAGE
        To use this module, you must:
        import the package:
            import sequenceAnalysis
        and call each of the functions:
            sequenceAnalysis.open_file(<<parameter>>)
FUNCTIONS
        open_file(<<file_name>>):
            return file_to_work
        create_fasta(<<file>>, <<sequences>>):
"""


def open_file(arguments):

    """
    This function takes a fasta file as input, opens it,
    reads all the lines, and stores them in a list.

    :param arguments: .fasta file
        Is the fasta file entered by command line
    :return file_to_work: list
        List with all the lines of the file
    """

    try:
        input_file = open(arguments)
    except IOError as ex:
        print("Sorry, couldn't find the file: " + ex.strerror)
        alternative_file = input("Please enter your file: ")
        input_file = open(alternative_file)

    file_to_work = input_file.readlines()
    return file_to_work
    input_file.close()


def create_fasta(file, sequences):

    """
    This function creates a .fasta file and stores
    any type of sequences in the same format

    :param file: .fasta file
        Is the path of the output FASTA file
        entered by command line
    :param sequences: list
        List with all the sequences
    """

    fasta_file = open(file, "x")
    for i, sequence in enumerate(sequences):

        # The head of each line is added to the fasta file with a '>' symbol, along with a new line.
        head = "> " + "seq_" + str(i + 1)
        fasta_file.write(head + "\n")

        # The DNA sequences are added to the fasta file.
        fasta_file.write(sequence.upper() + "\n")

    # The fasta file is closed
    fasta_file.close()

    # A message confirming the creation of the file
    print("The fasta file was created successfully! Check the output path provided.")

