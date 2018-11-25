results energies_and_pressures(state molA, state molB, double &Lx, double &Ly, double &beta)
{

results en_and_press = energies_and_forces_2(molA, molB, Lx, Ly, beta);
results testing_press = en_and_press;

double x_distance = molA.x-molB.x;
double y_distance = molA.y-molB.y;

state mol1 = molA;
state mol2 = molB;

double dr = 0.01;
mol1.x = PBC2D(Lx, molA.x + dr);
mol2.x =  PBC2D(Lx, molB.x - dr);
results energy_one = energies_and_forces_2(mol1, mol2, Lx, Ly,beta);

mol1.x = PBC2D(Lx, molA.x - dr);
mol2.x = PBC2D(Lx, molB.x + dr);
results energy_two = energies_and_forces_2(mol1, mol2, Lx, Ly,beta);

en_and_press.p_X_LJ = -(energy_one.energy_LJ - energy_two.energy_LJ)/(dr*4.0)*x_distance;
en_and_press.p_X_QQ = -(energy_one.energy_QQ - energy_two.energy_QQ)/(dr*4.0)*x_distance;

mol1 = molA;
mol2 = molB;

mol1.y = PBC2D(Ly, molA.y + dr);
mol2.y = PBC2D(Ly, molB.y - dr);
energy_one = energies_and_forces_2(mol1, mol2, Lx, Ly,beta);

mol1.y = PBC2D(Ly, molA.y - dr);
mol2.y = PBC2D(Ly, molB.y + dr);
energy_two = energies_and_forces_2(mol1, mol2, Lx, Ly,beta);

en_and_press.p_Y_LJ = -(energy_one.energy_LJ - energy_two.energy_LJ)/(dr*4.0)*y_distance;
en_and_press.p_Y_QQ = -(energy_one.energy_QQ - energy_two.energy_QQ)/(dr*4.0)*y_distance;

if (testing_press.p_X_LJ/en_and_press.p_X_LJ > 1.3 || testing_press.p_X_LJ/en_and_press.p_X_LJ < 0.7)
{
cout << "Coordinates_1:" << molA.x << " " << molA.y << " " << molA.phi << endl;
cout << "Coordinates_2:" << molB.x << " " << molB.y << " " << molB.phi << endl;
cout << "PRESS_ONE:" << testing_press.p_X_LJ << " " << testing_press.p_X_QQ << " " << testing_press.p_Y_LJ << " " << testing_press.p_Y_QQ << endl;
cout << "PRESS_TWO:" << en_and_press.p_X_LJ << " " << en_and_press.p_X_QQ << " " << en_and_press.p_Y_LJ << " " << en_and_press.p_Y_QQ << endl;
}
return en_and_press;
}
