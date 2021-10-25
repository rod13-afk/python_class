"""
NAME
       arraysEstructurados.py
VERSION
        [1.0]
AUTHOR
        Rodrigo Daniel Hernandez Barrera <<rodrigoh@lcg.unam.mx>>
DESCRIPTION
        Este script imprime los arrays estructurados de los arrays vistos en
        clase, con el uso de la libreria numpy.
CATEGORY
        Numpy.
INPUT
        Dos arrays
OUTPUT
        Dos arrays estructurados
EXAMPLES
        Input:
            produccion = np.array([[5, 3], [11, 7], [4, 9], [2, 6]])
            costos = np.array([3.5, 5, 7, 4.3])
        Output:
            [('Gen_1',  5, 3) ('Gen_2', 11, 7) ('Gen_3',  4, 9) ('Gen_4',  2, 6)]
            [('Gen_1', 3.5) ('Gen_2', 5. ) ('Gen_3', 7. ) ('Gen_4', 4.3)]

GITHUB
        https://github.com/rod13-afk/python_class/blob/master/Tareas/arraysEstructurados.py

"""


import numpy as np

# Arrays
produccion = np.array([[5, 3], [11, 7], [4, 9], [2, 6]])
print('Arrays \n', produccion)
costos = np.array([3.5, 5, 7, 4.3])
print(costos, '\n')

# Arrays estructurados
produccion_estruc = np.array([('Gen_1', 5, 3), ('Gen_2', 11, 7), ('Gen_3', 4, 9), ('Gen_4', 2, 6)],
                                   dtype=[('nombre de gen', (np.str_, 10)), ('30 Celsius', np.int32), ('35 Celsius', np.int32)])

costos_estruc = np.array([('Gen_1', 3.5), ('Gen_2', 5), ('Gen_3', 7), ('Gen_4', 4.3)],
                         dtype=[('nombre de gen', (np.str_, 10)), ('Costo de induccion', np.float64)])

print('Arrays estructurados \n', produccion_estruc)
print(costos_estruc)
