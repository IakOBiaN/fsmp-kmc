using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{

 double energy;

 for(int i = 0; i < nPart; i++){coordinates[i].energy = 0;}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
       energy = Inter_potential(coordinates[molA], coordinates[molB], Lx, Ly, beta);
       coordinates[molA].energy += energy;
       coordinates[molB].energy += energy;
      }
    }
}
