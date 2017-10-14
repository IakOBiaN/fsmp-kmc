using namespace std;

void initConfigHerringbone (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{
double ratio_x_to_y = 26.826/21.34;//26.826/21.34;
double area = nPart*1000000/(state_dens*pow(sigma,2)*N_a);
Lx = sqrt(area*ratio_x_to_y);
Ly = Lx/ratio_x_to_y;

double a = Lx/sqrt(nPart);                       // half of x-unit vector of herringbone unit cell in sigma
double b = Ly/sqrt(nPart);          // half of y-unit vector of herringbone unit cell in sigma
//double gamma = 88.5*(3.141592653589/180.0);    // angle between unit vectors of the herringbone unit cell

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
                          coordinates[molecule].phi = -45.0*(3.141592653589/180.0);
                          molecule++;
                         }
            }
            else
            {
            if((j%2)!=0){
                        coordinates[molecule].x = PBC2D(Lx, j*a/2.0);
                        coordinates[molecule].y = PBC2D(Ly, i*b);
                        coordinates[molecule].phi = 45.0*(3.141592653589/180.0);
                        molecule++;
                       }
            }
        }
    }

 nPart = molecule;
 density = nPart/(Lx*Ly)/pow(sigma,2)/N_a*1000000;

 cout << "Herringbone structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mikro mol/m2" << "\t" << "Lx and Ly: " << Lx << " and " << Ly << endl;
}
