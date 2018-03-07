using namespace std;

void initConfigRandomTMA (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &state_dens)
{
double ratio_x_to_y = 26.826/21.34; //good for 20 K and p = 10.5

if (abs(state_dens-10.1)<0.01) {ratio_x_to_y = 1/0.76973;}
if (abs(state_dens-10.2)<0.01) {ratio_x_to_y = 1/0.77238;}
if (abs(state_dens-10.3)<0.01) {ratio_x_to_y = 1/0.77565;}
if (abs(state_dens-10.4)<0.01) {ratio_x_to_y = 1/0.77783;}
if (abs(state_dens-10.5)<0.01) {ratio_x_to_y = 1/0.77832;}
if (abs(state_dens-10.6)<0.01) {ratio_x_to_y = 1/0.77997;}
if (abs(state_dens-10.7)<0.01) {ratio_x_to_y = 1/0.7828;}
if (abs(state_dens-10.8)<0.01) {ratio_x_to_y = 1/0.78274;}
if (abs(state_dens-10.9)<0.01) {ratio_x_to_y = 1/0.78328;}
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
                          coordinates[molecule].sin_phi = sin(coordinates[molecule].phi);
                          coordinates[molecule].cos_phi = cos(coordinates[molecule].phi);
                          molecule++;
                         }
            }
            else
            {
            if((j%2)!=0){
                        coordinates[molecule].x = PBC2D(Lx, j*a/2.0);
                        coordinates[molecule].y = PBC2D(Ly, i*b);
                        coordinates[molecule].phi = 45.0*(3.141592653589/180.0);
                        coordinates[molecule].sin_phi = sin(coordinates[molecule].phi);
                        coordinates[molecule].cos_phi = cos(coordinates[molecule].phi);
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
