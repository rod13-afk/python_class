"""
NAME
       descripciones&Articulos.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Este script contiene dos funciones. La primera de ellas imprime la descripcion
        de un campo de la base de datos 'protein'. La segunda funcion crea un archivo de
        texto que contiene los IDs de los articulos que concordaron con el termino de busqueda.
CATEGORY
        NCBI data bases.
INPUT
        id_articulos tiene como input:
             El path donde se guardara el output
             El nombre del autor que se va a buscar
             la primera palabra que se va a buscar
             la segunda palabra que se va a buscar
OUTPUT
        Output de la funcion 'descripciones':
            Imprime el nombre del field y su descripcion
        Output de la funcion 'id_articulos':
            Archivo de texto en el path especificado
EXAMPLES

        Funcion 'descripciones':
            Input:
                campos = [["FieldList", "ECNO"], ["LinkList", "protein_protein_small_genome"]]
            Output:
                ECNO ---> EC number for enzyme or CAS registry number
                protein_protein_small_genome ---> All proteins from this genome

        Funcion 'id_articulos':
            Input:
                Ingrese el path donde se guardara el output: ../output/articulos.txt
                Introduzca el nombre del autor: Giampaolo
                Introduzca la primera palabra: cancer
                Introduzca la segunda palabra: heterogeneous
            Output (dentro de archivo):
                Resultados:
                ['27343631', '27105536', '23667824', '23411686', '19708094', '19081479', '18848786', '18068853', '16110750']

GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/descripciones&Articulos.py

"""

from Bio import Entrez

# Tarea 1

# Correo
Entrez.email = "rodrigoh@lcg.unam.mx"

# Handle con einfo
handle = Entrez.einfo(db="protein")  # indicar db de interes, data base
record = Entrez.read(handle)


def descripciones(campo):

    """
    Esta funcion recibe el campo que se desea buscar y el
    nombre del field que se encuentra dentro de este. Imprime
    directamente la descripcion del field.

    :param campo:
    :return:
    """

    for elemento in record["DbInfo"]:
        if elemento == campo[0]:
            for field in record["DbInfo"][campo[0]]:
                if field["Name"] == campo[1]:
                    print(field["Name"] + " ---> " + field["Description"])


# Ejemplo donde se guarda cada campo con su respectivo field. Un for llama a
# la funcion que imprime la descripcion, pasando como parametro cada uno de los
# campos.
campos = [["FieldList", "ECNO"], ["LinkList", "protein_protein_small_genome"]]
for campo in campos:
    descripciones(campo)

handle.close()


# Tarea 2

def id_articulos(path, name, palabras):
    """
    Esta funcion crea un archivo de texto que contiene los
    id de cada articulo que cumple con los criterios de busqueda.

    :param path: Es una cadena que guarda el path del archivo de salida
    :param name: Es una cadena que guarda el nombre del autor que se va abuscar
    :param palabras: Es una lista que guarda las palabras que se van a buscar en el titulo
    :return: un archivo de texto que se guarda en el path indicado
    """

    outputFile = open(path, "x")

    termino = f'({name}[AUTH] AND ({palabras[0]}[Title] OR {palabras[1]}[Title]))'
    handle = Entrez.esearch(db="pubmed", term=termino)
    ids = Entrez.read(handle)

    outputFile.write("Resultados: \n" + str(ids["IdList"]))
    outputFile.close()


palabras = []
path = str(input("Ingrese el path donde se guardara el output: "))
name = str(input("Introduzca el nombre del autor: "))
palabras.append(input("Introduzca la primera palabra: "))
palabras.append(input("Introduzca la segunda palabra: "))

id_articulos(path, name, palabras)
