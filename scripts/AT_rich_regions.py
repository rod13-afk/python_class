"""
NAME
       AT_rich_regions.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Programa que toma una secuencia de DNA, evalua que esta secuencia contenga unicamente nucleotidos y solo asi
        imprime las secuencias ricas en AT, de lo contrario imprime las letras ambiguas encontradas y su posicion.
CATEGORY
        Genomic sequence.
INPUT
        Una secuencia de DNA
OUTPUT
        Regiones ricas en AT o la letra incorrecta con su posicion
EXAMPLES
        Input
            dna = CTGCATTATATCGTACGAAATTATACGCGCG
         Output
            Regiones ricas en AT:
            ATTATAT
            AAATTATA
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/scripts/AT_rich_regions.py
"""


import re


def search(dna):
    """
        Funcion que recibe una secuencia de DNA e imprime a pantalla las regiones ricas en AT que tienen mas de 5 A o T
    """

    # Expresion usada para identificar las letras no deseadas
    errores = re.finditer(r"[^ATGC]", dna)

    # Condicional que imprime las regiones ricas en AT si no encuentra letras distintas a los nucleotidos
    if not re.search(r"[^ATGC]", dna):
        AT_rich = re.finditer(r'(([AT]){5,})', dna)
        print("Regiones ricas en AT:")
        for intervalo in AT_rich:
            print(intervalo.group())
    else:

        # Imprime las letras ambiguas y su posicion
        for error in errores:
            base = error.group()
            pos = error.start()
            print("Se encontro", base, "en la posicion", pos)


dna = "CTGCATTATATCGTACGAAATTATACGCGCG"
search(dna)
