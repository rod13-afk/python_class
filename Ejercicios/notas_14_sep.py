# from Bio.Seq import Seq
from Bio import SeqIO

path = "../data/sample.fastq"
# Guardar ids de records bajo umbral
mala_calidad = []
umbral = 40
for record in SeqIO.parse(path, "fastq"):
    promedio = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
    if (promedio < umbral):
        mala_calidad.append((promedio, record.id, record.seq))

print(len(mala_calidad))

from Bio import SeqIO
for gb_record in SeqIO.parse("../data/aichi.gb", "genbank"):
    print('ID', gb_record.id)
    print('Secuencia', str(gb_record.seq)[0:30],'...')
    print('Longitud', len(gb_record))


path_2 = "../data/virus.gb"
for gb_record in SeqIO.parse(path_2, "genbank"):
    print("ID", str(gb_record.id))
print(gb_record.annotations)

version = gb_record.annotations["sequence_version"]
print(gb_record.annotations['organism'])

print(gb_record.features)
print(gb_record.features[0].qualifiers['country'])
print(gb_record.features[0].qualifiers['isolation_source'])




