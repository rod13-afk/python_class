"""
NAME
       getAbstracts.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Este script contiene dos funciones. La primera crea un archivo de texto que contiene
        los IDs de los articulos que concordaron con el termino de busqueda. La segunda funcion
        lee el archivo de IDs antes creado y crea uno nuevo, en este nuevo archivo se guardan los
        abstracts de los articulos y tambien las citas de cada uno.
CATEGORY
        NCBI data bases.
INPUT
        id_articulos tiene como input:
             El path donde se guardara el output
             El nombre del autor que se va a buscar
             la primera palabra que se va a buscar
             la segunda palabra que se va a buscar

        get_abstract tiene como input
             El path donde se guardaran los abstracts

OUTPUT
        Output de la funcion 'id_articulos':
            Archivo de texto en el path especificado con IDs de articulos
        Output de la funcion 'get_abstract':
            Archivo de texto en el path especificado con abstracts e IDs de las citas de los articulos
EXAMPLES

        Funcion 'id_articulos':
            Input:
                Ingrese el path donde se guardara el output: ../output/articulos.txt
                Introduzca el nombre del autor: Giampaolo
                Introduzca la primera palabra: cancer
                Introduzca la segunda palabra: heterogeneous
            Output (dentro de archivo):
                Resultados:
                ['27343631', '27105536', '23667824', '23411686', '19708094', '19081479', '18848786', '18068853', '16110750']

        Funcion 'get_abstract'
            Input:
                Ingrese el path donde se guardaran los abstracts: ../output/abstracts.txt
            Output:
                <<Archivo con abstracts y los IDs de las citas de los articulos>>

GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/getAbstracts.py

"""

from Bio import Entrez


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


def get_abtract(file, path2):

    """
    Esta funcion recibe una variable 'file' que contiene
    los datos del archivo que almacena los ids, producto
    de la funcion anterior

    :param file: lista de lista que guarda las lineas del archivo
    :param path2: string que guarda el path para el nuevo archivo
    :return: crea un archivo txt con los abstract y las citas de cada uno
    """

    # Obtencion de los abstracts
    ids = str(file[1])
    out_handle = open(path2, "w")
    fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)

    # Obtencion de las citas de los articulos
    citas = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=ids))
    for i in enumerate(ids):

        pmc_ids = f'Las citas del articulo con ID {ids[i]} son: {[link["Id"] for link in citas[i]["LinkSetDb"][0]["Link"]]}'
        out_handle.write(pmc_ids)
        out_handle.write('\n')
        out_handle.close()


# Correo
Entrez.email = "rodrigoh@lcg.unam.mx"

palabras = []
path = str(input("Ingrese el path donde se guardara el output: "))
name = str(input("Introduzca el nombre del autor: "))
palabras.append(input("Introduzca la primera palabra: "))
palabras.append(input("Introduzca la segunda palabra: "))
id_articulos(path, name, palabras)

# Tarea 5

# Path para guardar los abstracts
path2 = str(input("Ingrese el path donde se guardaran los abstracts: "))

# Se abre el archivo anteriormente creado y se leen las lineas
file = open(path, 'r')
all_lines = file.readlines()
file.close()

# Se llama a la funcion con los parametros correspondientes.
get_abtract(all_lines, path2)










