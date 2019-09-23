using namespace std;

void initConfigHoneycombTMA_elongated (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{

double h_bond_dist = 9.975; // It requires r_min for h-bonding

double ratio_x_to_y = 3.0/sqrt(3); //We use symmetric box
double L_cell = 2*h_bond_dist*cos(30.0/180.0*PI); //size of the simple unit cell
double area = L_cell*2.0*L_cell*cos(30.0/180.0*PI); // Area of the surface in A^2 (one unit cell in x direction and two in y direction)
int steps = sqrt(nPart/4.0/(2*2))+1;
Lx = steps*L_cell*ratio_x_to_y; // Size of the simulation box in A^2
Ly = Lx/ratio_x_to_y;
Lx *= 8.0;

double x_step = L_cell*cos(30.0/180.0*PI);
double y_step = L_cell;
double y_dop_step = L_cell*sin(30.0/180.0*PI);

int molecule = 0; // Molecules counter

for(int i = 0; i < steps; i++)
    {
     for(int j = 0; j < steps*4.0; j++)
        {
          for (int k = 0; k < 2; k++)
          {
              //unit cell
              coordinates[molecule].x = PBC2D(Lx, Lx/2.0 + Lx/4.0 - j*2.0*x_step-k*x_step);
              coordinates[molecule].y = PBC2D(Ly, 5.555 + i*y_step + j*2.0*y_dop_step + k*y_dop_step);
              coordinates[molecule].phi = 90.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;

              coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x - h_bond_dist*sin(30.0/180.0*PI));
              coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y + h_bond_dist*cos(30.0/180.0*PI));
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;

        }
      }
    }

while(molecule > nPart)
  {
    int max_numb,min_numb;double max_coord = -1e20,min_coord = 1e20;
    for(int j = 0; j < molecule; j++)
    {
      if (coordinates[j].x > max_coord) {max_numb = j; max_coord = coordinates[j].x;}
      if (coordinates[j].x < min_coord) {min_numb = j; min_coord = coordinates[j].x;}
    }
    for(int k = max_numb; k < molecule-1; k++)
    {
      coordinates[k] = coordinates[k+1];
    }
    molecule--;
    if (max_numb < min_numb) {min_numb--;}  //because we removed one molecule
    if(molecule == nPart) {break;}
    for(int k = min_numb; k < molecule-1; k++)
    {
      coordinates[k] = coordinates[k+1];
    }
    molecule--;
  }

assert(nPart == molecule);

 density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2 NOT CORRECT!!!

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << "mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
