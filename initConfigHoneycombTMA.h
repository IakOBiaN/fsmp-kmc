using namespace std;

void initConfigHoneycombTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double state_Ly)
{
// collision diameter = 11.052 A
// Ly = 23.84 * 11.052 = 263.4797 A
// Lx = Lx*2/sqrt(3) = 27.53 = 27.53 * 11.052 =  304.24 A
// h_bond_dist = 10.97832 A
double h_bond_dist = 10.97832;

double ratio_x_to_y = 2.0/sqrt(3);
double x_uc = 2.0*h_bond_dist*cos(30.0/180.0*PI);
double y_uc = 2.0*h_bond_dist + 2.0*h_bond_dist*sin(30.0/180.0*PI);

int cells = nPart/4.0;
int number_in_y = sqrt(cells/2.0);
int number_in_x = cells/number_in_y;

Lx = number_in_x*x_uc;
Ly = number_in_y*y_uc;

int molecule = 0; // Molecules counter

for(int i = 0; i < number_in_x; i++)
    {
     for(int j = 0; j < number_in_y; j++)
        {
          //unit cell
          coordinates[molecule].x = 3.0 + i*x_uc;
          coordinates[molecule].y = 3.0 + j*y_uc;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule]);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule]);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x;
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule]);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
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
		coordinates[i].x *= Lx_correction;
		coordinates[i].y *= Ly_correction;
		charges_coordinates(coordinates[i]);
	}
Ly *= Ly_correction;
Lx *= Lx_correction;
}


density = (1.0e+26)*nPart/(Lx*Ly)/N_a;

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
