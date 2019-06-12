using namespace std;

void initConfigHexTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{

double h_bond_dist = 9.975; // It requires r_min for h-bonding

double ratio_x_to_y = (h_bond_dist+cos(60.0/180.0*PI)*h_bond_dist)/(sin(60.0/180.0*PI)*h_bond_dist); //We use symmetric box

double area = (1.0e+26)*nPart/(state_dens*N_a); // Area of the surface in A^2
Lx = sqrt(area*ratio_x_to_y); // Size of the simulation box in A^2
Ly = Lx/ratio_x_to_y;


// half of x-unit vector of herringbone unit cell in sigma
double a = cos(60.0/180.0*PI)*h_bond_dist*2.0+h_bond_dist*2.0;
// half of y-unit vector of herringbone unit cell in sigma
double b = sin(60.0/180.0*PI)*h_bond_dist;

int steps = sqrt(nPart);

int molecule = 0; // Molecules counter

double move = 0.0;

 for(int i = 0; i < steps; i++)
    {
     for(int j = 0; j < steps/2.0; j++)
        {
            if ((i%2)==0)
            {
              coordinates[molecule].x = PBC2D(Lx, move+j*a);
              coordinates[molecule].y = PBC2D(Ly, move+i*b);
              coordinates[molecule].phi = 0.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;

              coordinates[molecule].x = PBC2D(Lx, move+j*a+h_bond_dist);
              coordinates[molecule].y = PBC2D(Ly, move+i*b);
              coordinates[molecule].phi = 60.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;
            }
            else
            {
              coordinates[molecule].x = PBC2D(Lx, move+h_bond_dist+cos(60.0/180.0*PI)*h_bond_dist+j*a);
              coordinates[molecule].y = PBC2D(Ly, move+i*b);
              coordinates[molecule].phi = 0.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;

              coordinates[molecule].x = PBC2D(Lx, move+h_bond_dist+cos(60.0/180.0*PI)*h_bond_dist+j*a+h_bond_dist);
              coordinates[molecule].y = PBC2D(Ly, move+i*b);
              coordinates[molecule].phi = 60.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;
            }
        }
    }

 density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << "mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
