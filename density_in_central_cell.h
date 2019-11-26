using namespace std;

void density_in_central_cell (int &nPart, double &density, double &gas_density, double &centralPart, vector <state> &coordinates, double &Lx, double &Ly)
{
	double gas_max_coord = (Lx/4.0)*3.0, gas_min_coord = Lx/4.0;
	double mol_in_central_cell = 0, mol_in_gas = 0;

	for(int i = 0; i < nPart; i++)
	{
		if (coordinates[i].x < gas_min_coord || coordinates[i].x > gas_max_coord){mol_in_gas += 1.0;}
		mol_in_central_cell += coordinates[i].cent;
	}

 density = (1.0e+26)*mol_in_central_cell/((Lx/4.0)*Ly)/N_a; // Density in central cell
 gas_density = (1.0e+26)*mol_in_gas/((Lx/2.0)*Ly)/N_a; // Density in central cell
 centralPart = mol_in_central_cell;
}
