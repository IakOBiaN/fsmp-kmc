using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{
	results en_and_press;
	// Loop over all distinct particle pairs
	for(int mol = 0; mol < nPart; mol++)
	{
		coordinates[mol].en_and_pr = en_and_press;
	}
	for(int molA = 0; molA < (nPart - 1); molA++)
	{
		for(int molB = (molA + 1); molB < nPart; molB++)
			{
				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
				en_and_press = en_and_press/2.0;  //for molecules pair to value per molecule
				coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
				coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
			}
		coordinates[molA].ex_field_coeff = external_field(coordinates[molA].x, Lx);
		coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + coordinates[molA].ex_field_coeff;
		cout << "mol #" << molA << " ext: " << coordinates[molA].ex_field_coeff.energy << " u: " << coordinates[molA].en_and_pr.energy << endl;

	}
	coordinates[nPart-1].ex_field_coeff = external_field(coordinates[nPart-1].x, Lx);
	coordinates[nPart-1].en_and_pr = coordinates[nPart-1].en_and_pr + coordinates[nPart-1].ex_field_coeff;
	cout << "mol #" << nPart-1 << " ext: " << coordinates[nPart-1].ex_field_coeff.energy << " u: " << coordinates[nPart-1].en_and_pr.energy << endl;
}
