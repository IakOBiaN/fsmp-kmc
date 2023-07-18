using namespace std;

void initConfigHoneycombTMA_elongated_cell(int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double state_dens)
{

int cells = nPart/4.0/8.0;  //4 molecules in unit cell (amount)
int number_in_x = sqrt(cells) + 1;  //unit cells along x directions
int number_in_y = number_in_x;  //unit cells along y directions
nPart = number_in_x*number_in_y*32.0;
double state_Ly = sqrt((1.0e+26)*nPart/8.0/(state_dens*N_a/sqrt(3.0)));  //y-size of the cell from density
double y_uc = state_Ly/number_in_y;  //y-size of the unit cell
double x_uc = y_uc/sqrt(3.0);  //x-size of the unit cell
double h_bond_dist = y_uc/3.0;  //hydrogen bond distance related
number_in_x *= 8;

Lx = number_in_x*x_uc; //add vacuum slab at both sides
Ly = number_in_y*y_uc;  //y-size of the cell

int molecule = 0; // Molecules counter
for(int i = 0; i < number_in_x; i++)
    {
     for(int j = 0; j < number_in_y; j++)
        {
          //unit cell
         coordinates[molecule].x = i*x_uc;
          coordinates[molecule].y = j*y_uc;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x;
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
          molecule++;
      }
    }

    for (int i = 0; i < molecule; i++)
		{
			double abs_x;
			if (coordinates[i].x > Lx/2.0){abs_x = coordinates[i].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[i].x;}
			double ksi = 32.0*abs_x/Lx;
			if (ksi > 8.0)
			{
				coordinates[i] = coordinates[molecule-1];
				molecule--;
				i--;
			}
		}
		nPart = molecule;

/*
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
*/

 cout << endl << "Honeycomb TMA Structure in central cell: " << endl;
 cout << "N: " << molecule << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
