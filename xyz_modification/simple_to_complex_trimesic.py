from math import atan, degrees, sin, cos, radians

def print_molecule(f, xy, ang):
    bond_ring = 1.438
    bond_1 = 1.348
    bond_2 = 0.980
    bond_2 = 1.254
    molecule = []
    #ring
    x = xy[0] + bond_ring
    y = xy[1]
    molecule.append([x, y, "C"])
    x = x - sin(radians(30)) * bond_ring
    y = y - cos(radians(30)) * bond_ring
    molecule.append([x, y, "C"])
    x = x - bond_ring
    y = y
    molecule.append([x, y, "C"])
    x = x - sin(radians(30)) * bond_ring
    y = y + cos(radians(30)) * bond_ring
    molecule.append([x, y, "C"])
    x = x + sin(radians(30)) * bond_ring
    y = y + cos(radians(30)) * bond_ring
    molecule.append([x, y, "C"])
    x = x + bond_ring
    y = y
    molecule.append([x, y, "C"])
    
    #carboxyl
    x = xy[0] + 2 * bond_ring
    y = xy[1]
    molecule.append([x, y, "C"])
    x = x + sin(radians(56.3)) * bond_ring
    y = y + cos(radians(56.3)) * bond_ring
    molecule.append([x, y, "C"])
    for element in molecule:
        x = element[0]#xy[0] + (element[0] - xy[0]) * cos(radians(ang)) - (element[1] - xy[1]) * sin(radians(ang))
        y = element[1]#xy[1] + (element[0] - xy[0]) * sin(radians(ang)) + (element[1] - xy[1]) * cos(radians(ang))
        print(element[2], x, y, 0, flush = True, file = f)
    return 0

xyz_file = open("statistics.xyz", "r")
atoms = int(xyz_file.readline())
parameters_line = xyz_file.readline()

new_xyz_file = open("statistics_new.xyz", "w")
new_atoms = int(atoms / 4 * 8)
print(new_atoms, flush = True, file = new_xyz_file)
print(parameters_line, end = "", flush = True, file = new_xyz_file)

counter = 0
general_counter = 0

for line in xyz_file:
    counter += 1
    general_counter += 1
    if counter == 1:
        xy1 = list(map(float, line.split()[1:]))
    if counter == 4:
        counter = 0
        xy2 = list(map(float, line.split()[1:]))
        a = xy2[1] - xy1[1]
        b = xy2[0] - xy1[0]
        if b == 0:
            b = 1e-13
        angle = degrees(atan(a / b))
        if b < 0 and a > 0:
            angle = 180 + angle
        elif b < 0 and a < 0:
            angle = -180 + angle
        
        print_molecule(new_xyz_file, xy2, angle)
        
    if general_counter == atoms:
        exit()
        xyz_file.readline()
        xyz_file.readline()
        print(new_atoms, flush = True, file = new_xyz_file)
        print(parameters_line, end = "", flush = True, file = new_xyz_file)
        general_counter = 0
    
xyz_file.close()
new_xyz_file.close()