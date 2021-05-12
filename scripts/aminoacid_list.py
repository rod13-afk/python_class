'''
NAME
       Program that calculates the percentage of a list of amino acids.
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        This program takes a sequence of amino acids from a protein and a list of amino acids, looks for these in the
        sequence and returns the percentage, it also tests the robustness of the code, all done with functions.
CATEGORY
        Protein sequence.
INPUT
        A protein sequence or nothing if the user uses the example
OUTPUT
        Returns the percentage of the amino acid list in the sequence as output.
EXAMPLES
        Input
            MSRSLLLRFLLFLLLLPPLP
         Output
            The percentage of hydrophilic amino acids in the example sequence is: 65.0
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/scripts/aminoacid_list.py

'''


# This function tests the robustness of the code
def testing(protein):
    assert get_aa_percentage(example_protein, ["M"]) == 5
    assert get_aa_percentage(example_protein, ['M', 'L']) == 55
    assert get_aa_percentage(example_protein, ['F', 'S', 'L']) == 70
    assert get_aa_percentage(example_protein, hydrophilic_aa) == 65


# This function calculates the percentage of the amino acid list in the sequence
def get_aa_percentage(protein, aa_list):
    percentage = 0
    for aa in aa_list:
        percentage += protein.count(aa) * 100 / len(protein)
    return percentage


example_protein = "MSRSLLLRFLLFLLLLPPLP"
hydrophilic_aa = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']

# The program asks if the user wants to use a new protein or the provided example.
print("Do you want to use the protein from the example in class or do you want to enter a new protein?:  ")
print("To use the example enter 0")
print("To enter another protein enter 1")
decision = int(input("Decision: "))

# This 'if' evaluates the user's decision and calls the functions with the correct parameters
if decision == 0:
    print("The percentage of hydrophilic amino acids in the example sequence is: " + str(
        get_aa_percentage(example_protein, hydrophilic_aa)))
    testing(example_protein)

elif decision == 1:
    print("Enter the new protein: ")
    new_protein = input("Protein: ")
    print("The percentage of hydrophilic amino acids in the new sequence is: " + str(
        get_aa_percentage(new_protein.upper(), hydrophilic_aa)))
    testing(new_protein)
