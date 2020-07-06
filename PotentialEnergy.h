using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{
	results en_and_press;
	EN_AND_PR_counter = en_and_press;
	// Loop over all distinct particle pairs
	for(int molA = 0; molA < (nPart - 1); molA++)
	{
		for(int molB = (molA + 1); molB < nPart; molB++)
			{
				en_and_press = energies_and_forces_exact(coordinates[molA], coordinates[molB], Lx, Ly, beta);
//				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
				EN_AND_PR_counter = EN_AND_PR_counter + en_and_press;
			}
	}
}
