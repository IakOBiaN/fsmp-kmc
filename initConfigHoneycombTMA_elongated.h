using namespace std;

void initConfigHoneycombTMA_elongated (int &nPart, double &density, double &gas_density, double &centralPart, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{

double h_bond_dist = 9.74; // It requires r_min for h-bonding

double x_uc = 2.0*h_bond_dist*cos(30.0/180.0*PI);
double y_uc = 2.0*h_bond_dist + 2.0*h_bond_dist*sin(30.0/180.0*PI);


double free_space = 0.45;

int cells = nPart/4.0;
int number_in_y = 3;
int number_in_x = cells/number_in_y;

double pre_Lx = number_in_x*x_uc;

Lx = pre_Lx*(1.0 + free_space*2.0);
Ly = number_in_y*y_uc;

double move_from_border = (Lx - pre_Lx)/2.0;

int molecule = 0; // Molecules counter

for(int i = 0; i < number_in_x; i++)
    {
     for(int j = 0; j < number_in_y; j++)
        {
          //unit cell
          coordinates[molecule].x = 3.0 + move_from_border + i*x_uc;
          coordinates[molecule].y = 3.0 + j*y_uc;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].cent = cent_potential(coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].cent = cent_potential(coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x;
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].cent = cent_potential(coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].cent = cent_potential(coordinates[molecule].x, Lx);
          molecule++;
      }
    }
    nPart = molecule;
    cout << "MOLECULES:" << nPart << endl;
/*while(molecule > nPart)
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
  }*/


density_in_central_cell (nPart, density, gas_density, centralPart, coordinates, Lx, Ly);

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "gas density: " << gas_density << " mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
