using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{
	results en_and_press;
	// Loop over all distinct particle pairs
	for(int mol = 0; mol < (nPart - 1); mol++)
	{
		coordinates[mol].en_and_pr = en_and_press;
	}
	for(int molA = 0; molA < (nPart - 1); molA++)
	{
		for(int molB = (molA + 1); molB < nPart; molB++)
			{
				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
				en_and_press.energy /= 2.0;  //for molecules pair (pressure should be full)
				coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
				coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
			}
	}
}
