using namespace std;

void initConfigHoneycombTMA_elongated (int &nPart, double &density, double &gas_density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{

double h_bond_dist = 9.975; // It requires r_min for h-bonding

double y_uc = 2.0*h_bond_dist*cos(30.0/180.0*PI);
double x_uc = 2.0*h_bond_dist + 2.0*h_bond_dist*sin(30.0/180.0*PI);
int num_cells = nPart/4.0; // Number of unit cells (4 is the number of molecules in unit cell)





double x_step = L_cell*cos(30.0/180.0*PI);
double y_step = L_cell;
double y_dop_step = L_cell*sin(30.0/180.0*PI);

int molecule = 0; // Molecules counter

for(int i = 0; i < steps; i++)
    {
     for(int j = 0; j < steps; j++)
        {
          for (int k = 0; k < 2; k++)
          {
              //unit cell
              coordinates[molecule].x = j*2.0*x_step-k*x_step;
              coordinates[molecule].y = 5.555 + i*y_step + j*2.0*y_dop_step + k*y_dop_step;
              coordinates[molecule].phi = 0.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
              molecule++;

              coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*sin(30.0/180.0*PI);
              coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*cos(30.0/180.0*PI);
              coordinates[molecule].phi = 60.0;
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


density_in_central_cell (nPart, density, gas_density, coordinates, Lx, Ly);

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "gas density: " << gas_density << " mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
