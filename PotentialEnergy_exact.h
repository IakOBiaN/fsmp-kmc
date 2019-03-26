using namespace std;

void PotentialEnergy_exact(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{

 results en_and_press;

EN_AND_PR_counter = en_and_press;

if (rosenbluth) {for(int i = 0; i < nPart; i++) {coordinates[i].energy = 0;}}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
           en_and_press = energies_and_forces_exact(coordinates[molA], coordinates[molB], Lx, Ly,beta);
           if (rosenbluth)
           {
                coordinates[molA].energy += en_and_press.energy;
                coordinates[molB].energy += en_and_press.energy;
           }
           EN_AND_PR_counter = EN_AND_PR_counter + en_and_press;
      }
    }
}
