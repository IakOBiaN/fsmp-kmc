from math import atan, degrees, sin, cos, radians

#file_names
xyz_file_name = "lambda0_0.612372_T300_2.xyz"
molecule_file_name = "trimesic_acid.xyz"

#read molecule from xyz file
molecule = []
mol_file = open(molecule_file_name, "r")

atoms_in_mol = int(mol_file.readline())
mol_file.readline()

for line in mol_file:
    atom = line.split()
    atom = [atom[0]] + list(map(float, atom[1:]))
    molecule.append(atom)

mol_file.close()


#write molecule to new xyz file (function)
def print_molecule(f, xy, ang):
    for element in molecule:
        x = xy[0] + element[1] * cos(radians(ang)) - element[2] * sin(radians(ang))
        y = xy[1] + element[1] * sin(radians(ang)) + element[2] * cos(radians(ang))
        print(element[0], x, y, 0, flush = True, file = f)
    return 0

#read initial xyz file for convertation
xyz_file = open(xyz_file_name, "r")
#number of atoms in initial xyz file
atoms = int(xyz_file.readline())
#read parameters line for new file
parameters_line = xyz_file.readline()

#create new xyz file for modified version
new_xyz_file = open(xyz_file_name[:-4] + "_new.xyz", "w")

#number of atoms in new xyz file
new_atoms = int(atoms / 4 * atoms_in_mol)
print(new_atoms, flush = True, file = new_xyz_file)
print(parameters_line, end = "", flush = True, file = new_xyz_file)

counter = 0
general_counter = 0

for line in xyz_file:
    counter += 1
    general_counter += 1
    #remember coordinates of first atom for angle calculation
    if counter == 1:
        xy1 = list(map(float, line.split()[1:]))
    #remember coordinates of fourth atom for cener coordinates
    if counter == 4:
        counter = 0
        xy2 = list(map(float, line.split()[1:]))
        a =  xy1[1] - xy2[1]
        b = xy1[0] - xy2[0]
        if b == 0:
            b = 1e-13
        angle = degrees(atan(a / b))
        if b < 0 and a > 0:
            angle = 180 + angle
        elif b < 0 and a < 0:
            angle = -180 + angle
        print_molecule(new_xyz_file, xy2, angle)

    #new frame
    if general_counter == atoms:
        xyz_file.readline()
        xyz_file.readline()
        print(new_atoms, flush = True, file = new_xyz_file)
        print(parameters_line, end = "", flush = True, file = new_xyz_file)
        general_counter = 0
    
xyz_file.close()
new_xyz_file.close()