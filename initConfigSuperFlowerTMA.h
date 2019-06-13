using namespace std;

void initConfigSuperFlowerTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{

double h_bond_dist = 9.975; // It requires r_min for h-bonding
double h3_bond_dist = 9.759; // It requires r_min for h3-bonding

double ratio_x_to_y = 3.0/sqrt(3); //We use symmetric box
double L_cell = h3_bond_dist; //size of the simple unit cell
double area = L_cell*2.0*L_cell*cos(30.0/180.0*PI); // Area of the surface in A^2 (one unit cell in x direction and two in y direction)
int steps = sqrt(nPart/2.0);
Lx = steps*L_cell;//sqrt(area*ratio_x_to_y)*steps; // Size of the simulation box in A^2
Ly = Lx*ratio_x_to_y;

double x_step = L_cell;
double y_step = L_cell*cos(30.0/180.0*PI);
double x_dop_step = L_cell*sin(30.0/180.0*PI);

int molecule = 0; // Molecules counter

for(int i = 0; i < steps; i++)
    {
     for(int j = 0; j < steps; j++)
        {
          for (int k = 0; k < 2; k++)
          {
              //unit cell
              coordinates[molecule].x = PBC2D(Lx, 5.555 + i*x_step + j*2.0*x_dop_step + k*x_dop_step);
              coordinates[molecule].y = PBC2D(Ly, 5.555 - j*2.0*y_step-k*y_step);
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;
        }
      }
    }
nPart = molecule;
density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << "mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
