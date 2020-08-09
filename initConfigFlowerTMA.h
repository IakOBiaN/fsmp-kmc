using namespace std;

void initConfigFlowerTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double state_Ly)
{

double h_bond_dist = 10.9792; // It requires r_min for h-bonding
double h3_bond_dist = 11.1334; // It requires r_min for h3-bonding

double ratio_x_to_y = 2.0/sqrt(3); // We use symmetric box
double x_uc = 2*h_bond_dist*cos(30.0/180.0*PI) + h3_bond_dist;
double y_uc = x_uc/ratio_x_to_y;
double cell_move = y_uc*tan(30.0/180.0*PI);

int cells = nPart/6.0;
int number_in_y = sqrt(cells);
int number_in_x = cells/number_in_y;

Lx = number_in_x*x_uc;
Ly = number_in_y*y_uc;

int molecule = 0; // Molecules counter
double start_position = 50.0;

for(int i = 0; i < number_in_x; i++)
    {
     for(int j = 0; j < number_in_y; j++)
        {
              //unit cell
              coordinates[molecule].x = PBC2D(Lx, start_position + i*x_uc + j*cell_move);
              coordinates[molecule].y = PBC2D(Ly, start_position + j*y_uc);
              coordinates[molecule].phi = 90.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
              molecule++;

              coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x + h3_bond_dist);
              coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y);
              coordinates[molecule].phi = 90.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
              molecule++;


							coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI));
              coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y - h_bond_dist*sin(30.0/180.0*PI));
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
              molecule++;

							coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x + h3_bond_dist*sin(30.0/180.0*PI));
              coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y - h3_bond_dist*cos(30.0/180.0*PI));
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
              molecule++;

							coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x - h3_bond_dist);
							coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y);
							coordinates[molecule].phi = 30.0;
							coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
							molecule++;

							coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI));
							coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI));
							coordinates[molecule].phi = 90.0;
							coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							charges_coordinates(coordinates[molecule]);
							molecule++;

      }
    }
nPart = molecule;

if (state_Ly > 0)
{
	double Ly_correction = state_Ly/Ly;
	double Lx_correction = state_Ly*ratio_x_to_y/Lx;
	for (int i = 0; i < nPart; i++)
	{
		coordinates[i].x *= Ly_correction;
		coordinates[i].y *= Lx_correction;
		charges_coordinates(coordinates[i]);
	}
	Ly *= Ly_correction;
	Lx *= Lx_correction;
}

density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2

 cout << "Flower-1 TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << "mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
