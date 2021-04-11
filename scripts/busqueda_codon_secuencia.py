'''
NAME
       Busqueda del Codón inicial y secuencia transcrita
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernández Barrera
DESCRIPTION
        Find the start codon and the transcribed sequence
CATEGORY
        Genomic sequence
INPUT
        Read a DNA sequence entered by the user
OUTPUT
        Returns as output the start and end positions of the sequence
        to be transcribed and its nucleotide sequence
EXAMPLES
    Input
     dna = 'AAGGTACGTCGCGCGTTATTAGCCTAAT'
    Output
     El codon AUG empieza en la posicion 4 y termina en 19, tomando en cuenta el 0 como
     la posicion del primer nucleotido.
     Fragmento de RNA que es transcrito representado en DNA es: TACGTCGCGCGTTAT
     Fragmento de RNA que es transcrito representado en RNA es: UACGUCGCGCGUUAU
     
'''


print('Introduzca la secuencia de DNA de interes:')
dna = input()  # La secuencia input se almacena en la variable dna


codon_inicio = 'TAC'
codon_termino = 'ATT'


'''El metodo str.find() devuelve el índice más bajo en el que se encuentre
 el codon de inicio, aplicado tambien para encontrar la posicion del codon de termino'''
inicio = dna.find(codon_inicio)
final = dna.find(codon_termino)


'''Se corta la secuencia para obtener la secuencia transcrita y se suma 2 para
incluir el codon de paro completo en el output'''
transcrito = dna[inicio:final + 2]


print('El codon AUG empieza en la posicion ' + str(inicio) + ' y termina en ' + str(final + 2) + ', tomando en cuenta el 0 como la posicion del primer nucleotido.')
print('Fragmento de RNA que es transcrito representado en DNA es: ' + transcrito)
print('Fragmento de RNA que es transcrito representado en RNA es: ' + transcrito.replace('T', 'U'))
