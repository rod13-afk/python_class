"""
NAME
       funcionPDB.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Esta funcion recibe un archivo pdb, el nombre (id) de una cadena proteica
        y el residuo que desea buscarse en la cadena. Regresa una lista que contiene
        a los residuos encontrados y el id de cada uno.
CATEGORY
        Secuencia proteica.
INPUT
        Archivo pdb, el caracter que identifica a la cadena y las tres letras
        que representan al aminoacido que desea buscarse.
OUTPUT
        Lista con el nombre de los residuos y el id de cada uno
EXAMPLES
        Input
             residuos_chain_A = residuos('../data/1kcw.pdb', 'A', 'LYS')
        Output
             [['LYS', 1], ['LYS', 3], ['LYS', 23], ['LYS', 24], ['LYS', 49], ['LYS', 50], ... ]

        *El input se modifica en el main code que esta al final del codigo*

GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/funcionPDB.py

"""

# Librerias
from Bio import PDB


def residuos(path, c_name, r_name):

    """
    Esta funcion recibe un archivo pdb, el nombre de una cadena proteica y el
    nombre de un residuo y regresa los residuos encontrados con su nombre y su id.

    :param path: str que guarda el path del archivo.
    :param c_name: str que es la letra, nombre de la cadena.
    :param r_name: str las tres letras representativas del aminoacido.
    :return: list de listas con el nombre del residuo y su id.
    """

    # Lista vacia que guardara los residuos
    residuos = []

    # Se crea el parser
    parser = PDB.PDBParser(QUIET=True)

    obj_struct = parser.get_structure('protein', path)

    # Para cada modulo se va a acceder a cada cadena, si se encuentra el id (nombre) dado para la
    # cadena se recorren los residuos de la misma, y se van a guardar si el residuo es el especificado
    # en los parametros en la llamada de la funcion. Se guardan con append en la lista vacia
    for modelo in obj_struct:
        for chain in modelo:
            if chain.id == c_name:
                for residuo in chain:
                    if residuo.get_resname() == r_name:
                        residuos.append([residuo.get_resname(), residuo.get_id()[1]])

    return residuos


# Main code
residuos_chain_A = residuos('../data/1kcw.pdb', 'A', 'LYS')
print(residuos_chain_A)
