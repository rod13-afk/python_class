'''
NAME
       Counting DNA Nucleotides
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera
GitHub
        https://github.com/rod13-afk/python_class/blob/master/scripts/contenido_ATCG.py
DESCRIPTION
        Count the number of the four nucleotides in a DNA sequence
CATEGORY
        Genomic sequence
INPUT
        Read a DNA sequence entered by the user
OUTPUT
        Returns as output the number of the four nucleotides
EXAMPLES
    Input
     dna = 'AAGGAUGTCGCGCGTTATTAGCCTAA'
    Output
     Results:

     Adenine = 7
     Cytosine = 5
     Guanine = 7
     Thymine = 6

'''


print('Enter a DNA sequence\n')
sequence = input()

number_of_A = sequence.count('A')
number_of_C = sequence.count('C')
number_of_G = sequence.count('G')
number_of_T = sequence.count('T')

print('\nResults:\n\n' + 'Adenine = ' + str(number_of_A) + '\nCytosine = ' + str(number_of_C) + '\nGuanine = ' + str(number_of_G) + '\nThymine = ' + str(number_of_T))
