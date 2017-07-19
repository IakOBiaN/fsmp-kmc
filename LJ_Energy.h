using namespace std;

double LJ_Energy(int &nPart, double &sigma, double &Lx, double &Ly, double &Rc, vector <state> &coordinates)
{

 double energy = 0;
 double Rc2 = Rc * Rc;

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
       if((pow((coordinates[molB].x - coordinates[molA].x), 2) + pow((coordinates[molB].y - coordinates[molA].y), 2)) > Rc2) {continue;}
       energy += LJ_inter_mols(coordinates[molA], coordinates[molB], Lx, Ly);
      }
    }
 // For efficiency, we will multiply by 4 only after summing
 energy = 4.0 * energy;
 return energy;
}
