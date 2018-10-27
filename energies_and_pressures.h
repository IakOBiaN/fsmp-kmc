results energies_and_pressures(state molA, state molB, double &Lx, double &Ly, double &beta)
{

results en_and_press = energies_and_forces(molA, molB, Lx, Ly, beta);
double x_distance = molA.x-molB.x;
double y_distance = molA.y-molB.y;

state mol1 = molA;
state mol2 = molB;

double dr = 0.01;
mol1.x = molA.x + dr;
mol2.x = molB.x - dr;
results energy_one = energies_and_forces(mol1, mol2, Lx, Ly,beta);

mol1.x = molA.x - dr;
mol2.x = molB.x + dr;
results energy_two = energies_and_forces(mol1, mol2, Lx, Ly,beta);

en_and_press.p_X_LJ = (energy_one.energy_LJ - energy_two.energy_LJ)/(dr*4.0)*x_distance;
en_and_press.p_X_QQ = (energy_one.energy_QQ - energy_two.energy_QQ)/(dr*4.0)*x_distance;
if (molA.x > molB.x)
{
  en_and_press.p_X_LJ *= -1.0;
  en_and_press.p_X_QQ *= -1.0;
}

mol1 = molA;
mol2 = molB;

mol1.y = molA.y + dr;
mol2.y = molB.y - dr;
energy_one = energies_and_forces(mol1, mol2, Lx, Ly,beta);

mol1.y = molA.y - dr;
mol2.y = molB.y + dr;
energy_two = energies_and_forces(mol1, mol2, Lx, Ly,beta);

en_and_press.p_Y_LJ = (energy_one.energy_LJ - energy_two.energy_LJ)/(dr*4.0)*y_distance;
en_and_press.p_Y_QQ = (energy_one.energy_QQ - energy_two.energy_QQ)/(dr*4.0)*y_distance;

if (molA.y > molB.y)
{
  en_and_press.p_Y_LJ *= -1.0;
  en_and_press.p_Y_QQ *= -1.0;
}

return en_and_press;
}
