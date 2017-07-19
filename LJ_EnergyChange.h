using namespace std;

double LJ_EnergyChange (int &nPart, double &sigma, double &Lx, double &Ly, double &Rc, vector <state> &coordinates,
                        state &trial_mol, int &trialPart)
{
 double deltaE = 0;
 double eNew = 0, eOld = 0;
 double Rc2 = Rc * Rc;

 // Loop over all particles and calculate
 // interaction with particle "trialPart"
 for(int otherPart = 0; otherPart < nPart; otherPart++)
    {
     // Make sure to skip self-interaction
     if(otherPart == trialPart){continue;}
     if((pow((coordinates[otherPart].x - trial_mol.x), 2) + pow((coordinates[otherPart].y - trial_mol.y), 2)) > Rc2) {continue;}

     // Claculate particle-particle distance for both
     // old nad new configurations
     eNew += LJ_inter_mols(trial_mol, coordinates[otherPart], Lx, Ly);
     eOld += LJ_inter_mols(coordinates[trialPart], coordinates[otherPart], Lx, Ly);
     deltaE = deltaE + (eNew - eOld);
    }
  // For efficiency, we will multiply by 4 only after summing
  deltaE = 4.0 * deltaE;
  return deltaE;
}
