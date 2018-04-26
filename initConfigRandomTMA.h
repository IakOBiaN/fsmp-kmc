using namespace std;

void initConfigRandomTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{
double ratio_x_to_y = 1.0; //now it is unknown

//if (abs(state_dens-10.1)<0.01) {ratio_x_to_y = 1/0.76973;}
double area = nPart/(state_dens*N_a);
Lx = sqrt(area*ratio_x_to_y);
Ly = Lx/ratio_x_to_y;

int steps = nPart;

int molecule = 0; // Molecules counter

 for(int i = 0; i < steps; i++)
    {
      coordinates[molecule].x = RanGen.Random()*Lx;
      coordinates[molecule].y = RanGen.Random()*Ly;
      coordinates[molecule].ang = RanGen.Random()*360.0;
      molecule++;
    }

 nPart = molecule;
 density = nPart/(Lx*Ly)/N_a;

 cout << "Random TPA structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly: " << Lx << " and " << Ly << endl;
}
