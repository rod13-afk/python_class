"""
NAME
       sequenceAnalysis
VERSION
        1.0
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Package containing two modules, dnaAnalysis and fastaFiles, the first module contains functions
        to work with DNA sequences, the second has functions to open and create fasta files.
CATEGORY
        Genomic sequences management.
USAGE
        Import the modules and the functions.
"""

print(f'Invoking __init__.py for {__name__}')

import sequenceAnalysis.dnaAnalysis
import sequenceAnalysis.fastaFiles

from sequenceAnalysis.fastaFiles import open_file, create_fasta
from sequenceAnalysis.dnaAnalysis import only_sequences, transcription, get_AT_content
