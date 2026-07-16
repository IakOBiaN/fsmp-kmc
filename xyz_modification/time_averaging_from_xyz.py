from math import atan, degrees, sin, cos, radians
from random import randint, random
import matplotlib.pyplot as plt
import numpy as np

def distance(atom1, atom2):
    return ((atom1[1]-atom2[1])**2 + (atom1[2]-atom2[2])**2 + (atom1[3]-atom2[3])**2)**0.5

#file_names
xyz_file_name = "statistics_T300_lambda0_0.612372.xyz"
molecule_file_name = "../models/trimesic_acid.xyz"

xy_matrix = []
random_point_per_bond = 20
scale = 10
frames_number = 1
frames_to_average = 200
shift = 100

rng = np.random.default_rng()
max_x = 0
max_y = 0
min_x = 1e10
lattice = [0, 0]

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

#create bonds list (atom_i_num, atom_j_num)
bonds = set()
for i, atom1 in enumerate(molecule):
    Bonds = []
    for j, atom2 in enumerate(molecule):
        if (i == j) or ((atom1[0] == "H") and (atom2[0] == "H")):
            continue
        if (distance(atom1, atom2) < 1.2) or (distance(atom1, atom2) < 1.61 and ((atom1[0] != "H") and (atom2[0] != "H"))):
            bonds.add(tuple(sorted((i, j))))

#if you want to check bonds in molecule
"""print(len(bonds))
#show molecule and calculated structure
plt.figure()
#show bonds
for bond in bonds:
    atom1 = molecule[bond[0]]
    atom2 = molecule[bond[1]]
    plt.plot([atom1[1], atom2[1]], [atom1[2], atom2[2]], "k-")

#show atoms
for i, atom in enumerate(molecule):
    color = "ko"
    if atom[0] == "C":
        color = "go"
    elif atom[0] == "H":
        color = "yo"
    elif atom[0] == "N":
        color = "go"
    elif atom[0] == "O":
        color = "ro"
    plt.plot(atom[1], atom[2], color)
    plt.text(atom[1] + 0.05, atom[2] + 0.1, str(i) + "_" + atom[0])
plt.show()
exit()"""

#write molecule to new xyz file (function)
def print_molecule(f, xy, ang):
    global xy_matrix
    molecule_paint = []
    global rng
    global max_x
    global max_y
    global min_x
    randoms = rng.random((random_point_per_bond,))
    for element in molecule:
        x = xy[0] + element[1] * cos(radians(ang)) - element[2] * sin(radians(ang))
        y = xy[1] + element[1] * sin(radians(ang)) + element[2] * cos(radians(ang))
        molecule_paint.append([x, y])
    for bond in bonds:
        x1 = molecule_paint[bond[0]][0]
        y1 = molecule_paint[bond[0]][1]
        x2 = molecule_paint[bond[1]][0]
        y2 = molecule_paint[bond[1]][1]
        a = x2 - x1
        b = y2 - y1
        for point in randoms:
            x = int(scale * (x1 + a * point))
            y = int(scale * (y1 + b * point))
            if max_x < x:
                max_x = x
            if min_x > x:
                min_x = x
            if max_y < y:
                max_y = y
            if y < 0:
                y += lattice[1]
            if y > lattice[1]:
                y -= lattice[1]
            xy_matrix[x][y] += 1
    return 0

#read initial xyz file for convertation
xyz_file = open(xyz_file_name, "r")
#number of atoms in initial xyz file
atoms = int(xyz_file.readline())
#read parameters line
parameters_line = xyz_file.readline()

lattice = list(map(float, parameters_line[9:-2].split()))
lattice = [int(scale * lattice[0]), int(scale * lattice[4])]


xy_matrix = np.zeros((10000, 2000))
#create dat file for time averaging
new_dat_file = open(xyz_file_name[:-4] + "_new.dat", "w")

counter = 0
general_counter = 0
frame = 0

for line in xyz_file:
    general_counter += 1
    if (frame >= frames_number) and (frame < (frames_number + frames_to_average)):
        if (((frame - frames_number) % (frames_to_average / 100)) == 0) and (general_counter == 1):
            print("Reading xyz:",((frame - frames_number) / (frames_to_average / 100)), "%")
        #remember coordinates of first atom for angle calculation
        counter += 1
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
            print_molecule(new_dat_file, xy2, angle)

    #new frame
    if general_counter == atoms:
        xyz_file.readline()
        xyz_file.readline()
        general_counter = 0
        frame += 1

max = xy_matrix.max()
for i in range(max_x + min_x):
    xy_matrix[i][0] = max
    xy_matrix[i][1] = max
    xy_matrix[i][lattice[1]] = max
    xy_matrix[i][lattice[1] + 1] = max

for i in range(max_y + shift):
    print("Writing dat file:", i / max_y * 100, " %")
    for j in range(max_x + min_x):
        print(xy_matrix[j][max_y + int(shift / 2) - 1 - i], end=" ", file = new_dat_file)
    print("", file = new_dat_file)

xyz_file.close()
new_dat_file.close()