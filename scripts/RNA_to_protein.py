"""
NAME
       Translating RNA into Protein
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Programa que traduce una secuencia de RNA en proteina
CATEGORY
        Genomic sequence.
INPUT
        Secuencia de RNA desde teclado
OUTPUT
        Secuencia aminoacidica
EXAMPLES
        Input
            AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
         Output
             ['M', 'A', 'M', 'A', 'P', 'R', 'T', 'E', 'I', 'N', 'S', 'T', 'R', 'I', 'N', 'G', '']
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/scripts/RNA_to_protein.py

"""


import re

gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T',
    'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AAC': 'N', 'AAT': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R',
    'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H',
    'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R',
    'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
    'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G',
    'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S',
    'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L',
    'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W'}

proteina = []
rna = input("Ingresa la secuencia de RNA: ")
rna = rna.upper()
dna = rna.replace('U', 'T')

# Se encuetra el codon de inicio y se guarda en una variable
cadena_codificante = re.search(r"ATG[ATGC]+", dna)

# Solo si existe el codon de inicio se traduce la secuencia a proteina
if cadena_codificante:

    # Con group extraemos el patron donde la secuencia comienza con el codon de inicio
    cadena_codificante = cadena_codificante.group()

    # Se separan los codones de 3 en 3 con un for que va desde la posicion 0 hasta la ultima representada por la
    # longitud de la secuencia
    codones = [cadena_codificante[i:i + 3] for i in range(0, len(cadena_codificante), 3)]

    # Este for itera sobre los codones para traducir cada uno a proteina
    for codon in codones:

        # Con un if se confirma que cada codon tenga longitud de 3, si no es asi es posible que la secuencia este
        # incompleta y se imprime el mensaje correspondiente
        if len(codon) != 3:
            print("El ultimo codon de tu secuencia no es un triplete, tu secuencia puede estar incompleta.")
        else:
            aa = gencode.get(codon)
            proteina.append(aa)

            # Si el amino acido es '', es decir un codon de paro, la traduccion se detiene
            if aa == '':
                print("Se encontro un codon de paro y se detuvo la traduccion del RNA a proteina.")
                break
            else:
                pass
    print("La proteina codificada por el RNA es:", proteina)
else:
    print("No existe codon de inicio.")
