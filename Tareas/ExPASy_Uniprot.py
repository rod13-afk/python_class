"""
NAME
       ExPASy_Uniprot.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Este script contiene una funcion. Esta funcion busca los terminos GO especificados en los archivos raw
        identificados con los ID que tambien se proporcionan como parametro de la funcion. Se crea un archivo, en el
        cual se guarda el ID y nombre de la proteina donde se encontro el ID, el termino GO encontrado y su definicion
        el organismo al que pertenece la proteina, su localizacion celular, un abstract de uno de los articulos
        mencionados en el archivo y por ultimo la documentacion de los prosite ID.
CATEGORY
        ExPASy & Uniprot.
INPUT
        Una lista con los terminos GO.
        Una lista con los ID de uniprot
        Una cadena de caracteres con el path del archivo

OUTPUT
        Archivo de texto que contiene la informacion
EXAMPLES
        Input:
            uniprot_IDs = ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC",
                 "POLG_YEFV1", "POLG_DEN1W",
                 "Q6W352_CVEN9", "D9SV67_CLOC7",
                 "A9KSF7_LACP7", "B8I7R6_RUMCH"]
            GO_Terms = ["GO:0046755", "GO:0046761",
                 "GO:0046760", "GO:0039702",
                 "GO:0046765", "GO:0046762"]

            path = "../output/GO_Uniprot.txt"
            GO_UNIPROT(GO_Terms, uniprot_IDs, path)

        Output (dentro de archivo):
                        .
                        .
                        .
                        .
                        .
                ID: A8R4D4_9BETC
                Name: Envelope small membrane protein
                GO: GO:0046760
                GO Definition: P:viral budding from Golgi membrane
                Organism: Equine coronavirus.
                Location:  Host Golgi apparatus membrane
                Abstract:
                [...]
                Documentation of PDOC51555:
                [...]
                        .
                        .
                        .
                        .
                        .

GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/ExPASy_Uniprot.py

"""


from Bio import ExPASy, SwissProt, Entrez
from Bio.ExPASy import Prosite, Prodoc


def GO_UNIPROT(GO_Terms, uniprot_IDs, path):

    """
    Funcion que busca los terminos GO dados en los archivos
    SwissProt de cada ID proporcionado, para posteriormente
    crear un archivo y guardar la informacion especificada
    en el header
    :param GO_Terms: una lista que contiene los terminos GO a buscar
    :param uniprot_IDs: una lista que contiene los ID de SwissProt
    :param path: una cadena de caracteres que tiene el path donde se guardara el archivo

    """

    # Se abre el archivo nuevo
    outputFile = open(path, "w")

    # Por cada ID de SwissProt se obtienen los archivos crudos
    for ID in uniprot_IDs:
        handle = ExPASy.get_sprot_raw(ID)
        record = SwissProt.read(handle)

        # Se accede al campo cross_references y, por cada elemento que se encuentra, se buscan los terminos GO
        for element in record.cross_references:
            for GO in GO_Terms:
                id_prosite = []
                pdocs = []

                # Si el termino GO se encuentra, se extrae y guarda la informacion
                if GO in element:

                    # Se obtiene el nombre
                    name = record.description.split("=")[1]
                    name = name.split("{")[0]
                    name = name.split(";")[0]

                    # Se guarda el ID y el nombre de la proteina
                    outputFile.write("ID: " + ID + "\n")
                    outputFile.write("Name: " + name + "\n")

                    # Se guarda el termino GO y su definicion
                    outputFile.write("GO: " + GO + "\n")
                    outputFile.write("GO Definition: " + str(element[2]) + "\n")

                    # Se guarda el organismo
                    outputFile.write("Organism: " + record.organism + "\n")

                    # Este for busca en el record.comments la localizacion celular y la guarda
                    for i in range(0, len(record.comments)):
                        if "SUBCELLULAR" in record.comments[i]:

                            location = record.comments[i].split("SUBCELLULAR LOCATION:")[1]
                            location = location.split("{")[0]
                            location = location.split(";")[0]
                            outputFile.write("Location: " + location + "\n")

                    # Se obtienen unas referencias de cada ID, esto debido a que solo se necesita 1 abstract. De una de
                    # las referencias se extrae el abstract y se guarda en el archivo
                    references = []
                    for reference in record.references[1:2]:
                        references.append(reference.references[0][1])

                        ids = str(references[0])
                        fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="text")
                        data = fetch_handle.read()
                        fetch_handle.close()
                        outputFile.write("Abstract: \n" + data + "\n")

                    # Se busca en el campo cross_references las referencias de prosite. En una lista se guardan los ID,
                    # estos ID buscan nuevamente con get_prosite_raw para obtener la documentacion de cada uno.
                    for reference in record.cross_references:
                        if 'PROSITE' in reference:
                            id_prosite.append(reference[1])

                            handle = ExPASy.get_prosite_raw(reference[1])
                            record = Prosite.read(handle)

                            pdocs.append(record.pdoc)

                    for pdoc in pdocs:
                        handle = ExPASy.get_prosite_raw(pdoc)
                        record = Prodoc.read(handle)

                        outputFile.write("Documentation of " + pdoc + ": \n" + record.text + "\n\n")

        outputFile.write("\n\n")
    outputFile.close()


# Codigo de ejemplo donde se especifica el codigo y se declaran las lista de los terminos GO y la de los ID de uniprot.
# Tambien se declara el path donde se va a crear el archivo y se llama a la funcion con sus parametros.
Entrez.email = "rodrigoh@lcg.unam.mx"

uniprot_IDs = ["A0A0K2RVI7_9BETC", "A8R4D4_9BETC",
                 "POLG_YEFV1", "POLG_DEN1W",
                 "Q6W352_CVEN9", "D9SV67_CLOC7",
                 "A9KSF7_LACP7", "B8I7R6_RUMCH"]
GO_Terms = ["GO:0046755", "GO:0046761",
              "GO:0046760", "GO:0039702",
              "GO:0046765", "GO:0046762"]

path = "../output/GO_Uniprot.txt"
GO_UNIPROT(GO_Terms, uniprot_IDs, path)
