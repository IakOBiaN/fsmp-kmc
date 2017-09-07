using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly,
                     double &Rc, double &Rc2, vector <state> &coordinates, const double &A, double &C_q, double &beta)
{

 double energy;

 for(int i = 0; i < nPart; i++){coordinates[i].energy = (-334.4/((1/beta)*36.4))*abs(sin(coordinates[i].tetta));}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
       energy = Inter_potential(coordinates[molA], coordinates[molB], Rc, Rc2, Lx, Ly, A, C_q, beta);
       coordinates[molA].energy += energy;
       coordinates[molB].energy += energy;
      }
    }
}
