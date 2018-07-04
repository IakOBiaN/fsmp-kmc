using namespace std;

void initConfigHerringbone (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{
double ratio_x_to_y = 1.2571; //good for 20 K and p = 10.5

double area = (1.0e+26)*nPart/(state_dens*N_a); // Area of the surface in A^2
Lx = sqrt(area*ratio_x_to_y); // Size of the simulation box in A^2
Ly = Lx/ratio_x_to_y;

double a = Lx/sqrt(nPart);          // half of x-unit vector of herringbone unit cell in sigma
double b = Ly/sqrt(nPart);          // half of y-unit vector of herringbone unit cell in sigma

int steps = sqrt(nPart);

int molecule = 0; // Molecules counter

 for(int i = 0; i < steps; i++)
    {
     for(int j = 0; j < steps*2; j++)
        {
            if ((i%2)==0)
            {
            if((j%2)==0){
                          coordinates[molecule].x = PBC2D(Lx, j*a/2.0);
                          coordinates[molecule].y = PBC2D(Ly, i*b);
                          coordinates[molecule].phi = -45.0;
                          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
                          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
                          molecule++;
                         }
            }
            else
            {
            if((j%2)!=0){
                        coordinates[molecule].x = PBC2D(Lx, j*a/2.0);
                        coordinates[molecule].y = PBC2D(Ly, i*b);
                        coordinates[molecule].phi = 45.0;
                        coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
                        coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
                        molecule++;
                       }
            }
        }
    }

 nPart = molecule;
 density = (1.0e+26)*nPart/(Lx*Ly)/N_a; // Density in mkMol/m2

 cout << "Herringbone structure of N2 on HOPG: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly: " << Lx << " and " << Ly << endl;
}
