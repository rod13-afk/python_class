"""
NAME
       dnaAnalysis
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This module works with DNA sequences, it has a function of preparing the sequences,
        discriminating between the sequences and their headers, in addition, it eliminates N's
        and other symbols. This function is called from other functions in the same module,
        these functions can transcribe the DNA sequences and get the AT content. In future
        versions this module may have more functions to complete the analysis of genomic sequences.
CATEGORY
        Genomic sequences.
USAGE
        To use this module, you must:
        import the package:
            import sequenceAnalysis
        and call each of the functions:
            sequenceAnalysis.transcription()
FUNCTIONS
        only_sequences(<<parameter>>)
            return sequences_to_work
        transcription(<<parameter>>)
            return RNA_sequences
        get_AT_content(<<parameter>>)
            print the sequence and its AT content
"""


def only_sequences(file):

    """
    First this function uses regular expressions to save only
    the sequences and not the fasta file headers, then the N's,
    quotes and hyphens are removed from the sequences.

    :param file: list
        This object contains all the lines of a fasta file
    :return: list
        Returns a list with the sequences ready to work with them
    """

    import re
    sequences = []
    for line in file:
        if re.search(r"[ATGC]", line):
            sequences.append(line)
        else:
            pass

    sequences_to_work = []
    for i, line in enumerate(sequences):

        sequence = str(line).replace('"', '')
        sequence = sequence.replace("N", '')
        sequence = sequence.replace('\n', '')
        sequences_to_work.append(sequence)

    return sequences_to_work


def transcription(file):

    """
    Function that transcribes a DNA sequence into RNA,
    first calls the function to prepare the sequences
    and then replaces the DNA nucleotides with their
    corresponding ones in RNA

    :param file: list
        This object contains all the lines of a fasta file
    :return: list
        Returns a list with the RNA sequences
    """

    RNA_sequences = []
    DNA_sequences = only_sequences(file)

    for sequence in DNA_sequences:
        sequence = sequence.replace("A", "U")
        sequence = sequence.replace("T", "a")
        sequence = sequence.replace("C", "g")
        sequence = sequence.replace("G", "c")
        RNA_sequences.append(sequence.upper())
    return RNA_sequences


def get_AT_content(file):

    """
    function that calculates the AT content by counting nucleotides A and T,
    adding them and dividing that sum by the length of the sequence

    :param file: list
        This object contains all the lines of a fasta file
    :return:
        prints the sequence number and its AT content
    """

    DNA_sequences = only_sequences(file)
    for i, sequence in enumerate(DNA_sequences):

        length = len(sequence)
        a_content = (sequence.upper()).count("A")
        t_content = (sequence.upper()).count("T")
        at_content = (a_content + t_content) / length

        print("seq_" + str(i + 1) + " has " + str(at_content) + " AT content" + "\n")


