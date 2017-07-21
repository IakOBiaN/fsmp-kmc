using namespace std;

double PotentialEnergy(int &nPart, double &sigma, double &Lx, double &Ly, double &Rc, vector <state> &coordinates, const double &A, double &C_q)
{

 double energy = 0;
 double Rc2 = Rc * Rc;

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
       if((pow((coordinates[molB].x - coordinates[molA].x), 2) + pow((coordinates[molB].y - coordinates[molA].y), 2)) > Rc2) {continue;}
       energy += Inter_potential(coordinates[molA], coordinates[molB], Lx, Ly, A, C_q);
      }
    }
 // For efficiency, we will multiply by 4 only after summing
 energy = 4.0 * energy;
 return energy;
}
