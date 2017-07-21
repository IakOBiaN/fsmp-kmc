using namespace std;

double initConfig (int &nPart, double &density, double &sigma, vector <state> &coordinates, double &beta, double &Rc, const double &A, double &C_q)
{
 double L; // Linear size of the box
 double Lx;
 double Ly;

 // Get the corresponding box size
 L = 2*nPart/density; // Lx*Ly
 Lx = sqrt(L*sqrt(3)/2);
 Ly = Lx*2/sqrt(3);

 double Rc2 = Rc*Rc;

 cout << "L: " << L << "\t" << "Lx: " << Lx << "\t" << "Ly: " << Ly << endl;

 int molecules = 0; // Molecules counter
 while(molecules < nPart)
  {
   state test_mol;
   test_mol.x = Lx * RanGen.Random();
   test_mol.y = Ly * RanGen.Random();
   test_mol.tetta = 2*3.141592653589 * RanGen.Random();
   test_mol.phi = 2*3.141592653589 * RanGen.Random();

   double U_LJ = 0;

   // Loop over all distinct particle pairs
   for(int mol = 0; mol < molecules; mol++)
     {
      if((pow((coordinates[mol].x - test_mol.x), 2) + pow((coordinates[mol].y - test_mol.y), 2)) > Rc2) {continue;}
      U_LJ += Inter_potential(test_mol, coordinates[mol], Lx, Ly, A, C_q);
     }
     U_LJ*=4;
     if(RanGen.Random() >= exp(-beta*U_LJ)){continue;}
     else
        {
            coordinates[molecules] = test_mol;
            molecules++;
        }
  }
 cout << "N: " << molecules << endl;
 return Lx;
}
