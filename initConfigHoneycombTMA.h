using namespace std;

void initConfigHoneycombTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly)
{
// collision diameter = 11.052 A
// Ly = 23.84 * 11.052 = 263.4797 A
// Lx = Lx*2/sqrt(3) = 27.53 = 27.53 * 11.052 =  304.2616 A
// h_bond_dist = 10.978 A
double h_bond_dist = 10.978;

double x_uc = 2.0*h_bond_dist*cos(30.0/180.0*PI); // 19.142 A
double y_uc = 2.0*h_bond_dist + 2.0*h_bond_dist*sin(30.0/180.0*PI); // 33.156 A

int cells = nPart/4.0;
int number_in_y = sqrt(cells/2.0);
int number_in_x = cells/number_in_y;

Lx = number_in_x*x_uc; // 306.272 if h_bond_dist = collision diameter
Ly = number_in_y*y_uc; // 265.248 if h_bond_dist = collision diameter


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
          charges_coordinates(coordinates[molecule], Lx, Ly);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule], Lx, Ly);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x;
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
          coordinates[molecule].phi = 30.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule], Lx, Ly);
          molecule++;

          coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
          coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
          coordinates[molecule].phi = 90.0;
          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
          charges_coordinates(coordinates[molecule], Lx, Ly);
          molecule++;
      }
    }
    nPart = molecule;
    cout << "MOLECULES:" << nPart << endl;
    density = (1.0e+26)*nPart/(Lx*Ly)/N_a;

 cout << "Honeycomb TMA Structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}
