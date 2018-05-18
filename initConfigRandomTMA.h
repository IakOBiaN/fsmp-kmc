using namespace std;

void initConfigRandomTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{
double ratio_x_to_y = 1.0; //We use symmetric box

//if (abs(state_dens-10.1)<0.01) {ratio_x_to_y = 1/0.76973;}
double area = (1.0e+26)*nPart/(state_dens*N_a); // Area of the surface in A^2
Lx = sqrt(area*ratio_x_to_y); // Size of the simulation box in A^2
Ly = Lx/ratio_x_to_y;

int steps = nPart;

int molecule = 0; // Molecules counter

 for(int i = 0; i < steps; i++)
    {
      coordinates[molecule].x = RanGen.Random()*Lx;
      coordinates[molecule].y = RanGen.Random()*Ly;
      coordinates[molecule].phi = RanGen.Random()*360.0;
      coordinates[molecule].sin_phi = sin(coordinates[molecule].phi);
      coordinates[molecule].cos_phi = cos(coordinates[molecule].phi);
      molecule++;
    }

 nPart = molecule;
 density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2

 cout << "Random TMA structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly: " << Lx << " and " << Ly << endl;
}
