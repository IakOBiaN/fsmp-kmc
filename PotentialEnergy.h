using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{
	results en_and_press;
	// Refresh the per-molecule field coefficients before any pair energy is
	// computed: the damping must be current (lambda0 changes with the
	// temperature point) and the pair energies below depend on it.
	for(int mol = 0; mol < nPart; mol++)
	{
		coordinates[mol].en_and_pr = en_and_press;
		coordinates[mol].damping_coeff = damping_field(coordinates[mol].x, Lx); // Lambda^1/2
		coordinates[mol].stat_weight = weights_for_central_cell (coordinates[mol].x, Lx);
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
		coordinates[molA].ex_field_coeff = external_field_and_mask(coordinates[molA].x, coordinates[molA].y, Lx);
		coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + coordinates[molA].ex_field_coeff;
//		cout << "mol #" << molA << " ext: " << coordinates[molA].ex_field_coeff.energy << " u: " << coordinates[molA].en_and_pr.energy << endl;
	}
	coordinates[nPart-1].ex_field_coeff = external_field_and_mask(coordinates[nPart-1].x, coordinates[nPart-1].y, Lx);
	coordinates[nPart-1].en_and_pr = coordinates[nPart-1].en_and_pr + coordinates[nPart-1].ex_field_coeff;
//	cout << "mol #" << nPart-1 << " ext: " << coordinates[nPart-1].ex_field_coeff.energy << " u: " << coordinates[nPart-1].en_and_pr.energy << endl;
}
