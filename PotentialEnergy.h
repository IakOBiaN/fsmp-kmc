using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{

results en_and_press_central_cell;
results en_and_press;

EN_AND_PR_counter = en_and_press;
EN_AND_PR_counter_central_cell = en_and_press_central_cell;

if (rosenbluth) {for(int i = 0; i < nPart; i++) {coordinates[i].energy = 0;}}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
	{
		for(int molB = (molA + 1); molB < nPart; molB++)
			{
				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta,true);
				//en_and_press_central_cell = en_and_press*coordinates[molA].cent*coordinates[molB].cent;

				EN_AND_PR_counter = EN_AND_PR_counter + en_and_press;
				//EN_AND_PR_counter_central_cell = EN_AND_PR_counter_central_cell + en_and_press_central_cell;
			}
	}
  EN_AND_PR_counter_central_cell = EN_AND_PR_counter;
}
