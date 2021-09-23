from Bio import PDB

parser = PDB.PDBParser(QUIET=True)
struc = parser.get_structure('prot_1fat', "../data/1fat.pdb")
print(struc.child_dict)
print(struc.child_list)

# Metadatos
print(struc.header.keys())
print(struc.header['structure_method'])
print(struc.header['resolution'])

# Ejercicio 1

parser = PDB.PDBParser(QUIET=True)

obj_struc = parser.get_structure("prot_1kcw", "../data/1kcw.pdb")
struct_method = obj_struc.header['structure_method']
struct_resolution = obj_struc.header['resolution']
print(struct_method, struct_resolution)

# Archivo de resonancia magnética nuclear
struc = parser.get_structure('RMN', "../data/1g03.pdb")
for model in struc:
    print(model)
    print(model.child_dict)

# Cadenas
for chain in model:
    print(chain)
    chain = model['A']
    print(chain)
    print(chain.child_dict)
    print(chain.child_list)

# Residuos
for residue in chain:
    print(residue)
    residuo = chain[45]
    print(residuo)
    print(residuo.get_id()[1])
    print(residuo.get_resname())
    residuos_int = []
    for residue in chain:
        if residuo.get_resname() == 'SER':
            residuos_int.append(residuo)
            print(len(residuos_int))

# Ejercicio 2
cisteinas = []
for modelo in obj_struc:
    for chain in modelo:
        if chain.id == 'A':
            for residuo in chain:
                if residuo.get_resname() == 'CYS':
                    cisteinas.append(residuo)
###########################################################################################

parser = PDB.PDBParser(QUIET=True)
obj_struc = parser.get_structure("prot_1kcw", "./files/clase_3/1kcw.pdb")
for key, valor in obj_struc.header.items():
    print(key, valor)
    struct_method = obj_struc.header['structure_method']
    struc_resolution = obj_struc.header['resolution']
    print(struct_method, struc_resolution)

    # Guardar las cisteínas
    cisteinas = []
    for modelo in obj_struc:
        for chain in modelo:
            if chain.id == 'A':
                for residuo in chain:
                    if residuo.get_resname() == 'CYS':
                        cisteinas.append(residuo)
                        print(len(cisteinas))
