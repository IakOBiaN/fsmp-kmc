using namespace std;

void initConfigSuperFlowerTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double state_Ly)
{

double h3_bond_dist = 11.1334; // It requires r_min for h3-bonding

double ratio_x_to_y = sqrt(3); //We use symmetric box
double x_uc = 2*h3_bond_dist*sin(60.0/180.0*PI);
double y_uc = x_uc/ratio_x_to_y;

int cells = nPart/2.0;
int number_in_y = sqrt(cells);
int number_in_x = cells/number_in_y;

Lx = number_in_x*x_uc;
Ly = number_in_y*y_uc;

int molecule = 0; // Molecules counter
double start_position = 3.0;

for(int i = 0; i < number_in_x; i++)
		{
     for(int j = 0; j < number_in_y; j++)
			{
              //unit cell
							coordinates[molecule].x = PBC2D(Lx, start_position + i*x_uc);
							coordinates[molecule].y = PBC2D(Ly, start_position + j*y_uc);
							coordinates[molecule].phi = 60.0;
							coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							molecule++;

							coordinates[molecule].x = PBC2D(Lx, coordinates[molecule-1].x + h3_bond_dist*sin(60.0/180.0*PI));
							coordinates[molecule].y = PBC2D(Ly, coordinates[molecule-1].y + h3_bond_dist*cos(60.0/180.0*PI));
							coordinates[molecule].phi = 60.0;
							coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
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
	}
	Ly *= Ly_correction;
	Lx *= Lx_correction;
}

density = (1.0e+26)*nPart/(Lx*Ly)/N_a;

cout << "Superflower TMA Structure: " << endl;
cout << "N: " << molecule << "\t" << "density: " << density << "mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
