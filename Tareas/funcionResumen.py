"""
NAME
       funcionResumen.py
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Esta funcion trabaja con archivos GenkBank, muestra un resumen
        de los metadatos y de la seccion 'feature', recibe una lista de
        genes que pueden ser buscados en el archivo, imprime el nombre
        de cada gen, los 15 primeros nucleotidos, asi como su respectiva
        secuencia de RNA y proteica.
CATEGORY
        GenBank

FUNCTIONS
        resumen(<path>, <genes>):
GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/funcionResumen.py
"""

from Bio import SeqIO

# Function
def resumen(path, genes):

    """
    Funcion que muestra la informacion almacenada
    en un archivo GenBank
    :param path: una cadena que especifica el path del archivo
    :param genes: una lista que guarda los genes que se van a buscar
    :return: no regresa ningun tipo de objeto, imprime todos los
    resultados dentro de la misma funcion.
    """



    for gb_record in SeqIO.parse(path, "genbank"):
        # Obtener e imprimir el organismo, buscando en los metadatos
        organismo = gb_record.annotations['organism']
        print('Organismo: ' + organismo)

        # Obtiene e imprime la fecha, buscando en los metadatos
        fecha = gb_record.annotations['date']
        print('Fecha: ' + fecha)

        # Obtiene e imprime el pais de la muestra, buscando en features
        pais = gb_record.features[0].qualifiers['country']
        print('Pais de la muestra: ' + str(pais))

        # Obtiene e imprime el numero del aislado, buscando en features
        aislado = gb_record.features[0].qualifiers['isolation_source']
        print('No. del aislado: ' + str(aislado) + '\n')

        # Uno de los for comienza con 1, porque el primer gen esta en la posicion 1 de features, lo hace de dos en dos
        # porque asi estan ordenados los genes
        for num in range(1, len(gb_record.features), 2):

            # El segundo for recorre los genes que se van a buscar
            for gene in genes:

                # El if es verdadero cuando encuentra al gen que estamos buscando, y procede a imprimir su nombre, la
                # secuencia de DNA, RNA y proteina de los primeros 15 nt
                if gb_record.features[num].qualifiers['gene'][0] == gene:

                    # Nombre de los genes
                    nombre = gb_record.features[num].qualifiers['gene']
                    print('Nombre del gen: ' + str(nombre))

                    # Los 15 nt del gene, su RNA y proteina
                    start = gb_record.features[num].location.nofuzzy_start
                    end = start + 15
                    DNA = gb_record.seq[start:end]
                    print(DNA)
                    RNA = DNA.transcribe()
                    print(RNA)
                    PROT = RNA.translate()
                    print(PROT, '\n')


# Main code, donde se especifica el path del archivo, los genes a buscar y se llama a la funcion.
path = "../data/virus.gb"
genes = ['L', 'N']
resumen(path, genes)
