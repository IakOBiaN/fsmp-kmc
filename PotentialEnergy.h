using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{

 results en_and_press;
 results en_and_press_2;

 double CC_max_coord = (Lx/16.0)*9.0, CC_min_coord = (Lx/16.0)*7.0;

EN_AND_PR_counter = en_and_press;

if (rosenbluth) {for(int i = 0; i < nPart; i++) {coordinates[i].energy = 0;}}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
           en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta,true);

					 if (rosenbluth)
           {
                coordinates[molA].energy += en_and_press.energy;
                coordinates[molB].energy += en_and_press.energy;
           }
           if ((coordinates[molA].x > CC_min_coord && coordinates[molA].x < CC_max_coord) || (coordinates[molB].x > CC_min_coord && coordinates[molB].x < CC_max_coord))
              {EN_AND_PR_counter = EN_AND_PR_counter + en_and_press;}
      }
    }
}
