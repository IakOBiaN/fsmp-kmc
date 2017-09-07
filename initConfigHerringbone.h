using namespace std;

void initConfigHerringbone (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly)
{

 double a = 4.04/3.318/2;                       // half of x-unit vector of herringbone unit cell in sigma
 double b = 7.00/3.318/2;                       // half of y-unit vector of herringbone unit cell in sigma
 double gamma = 88.5*(3.141592653589/180.0);    // angle between unit vectors of the herringbone unit cell

 // In herringbone structure there are 2 molecules per unit cell,
 // so we should include in the simulation box the following number of unit cells:

 int nUnitCells = sqrt(nPart/2.0);

 Lx = nUnitCells*a*2.0;
 Ly = nUnitCells*b*2.0;

 int stepsX = floor(Lx/a + 0.5);
 int stepsY = floor(Ly/b + 0.5);

 int molecule = 0; // Molecules counter

 for(int i = 0; i < stepsX; i++)
    {
     for(int j = 0; j < stepsY; j++)
        {
            if((i%2)==0)
            {
             if((j%2)==0){
                          coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                          coordinates[molecule].y = PBC2D(Ly, i*b);
                          coordinates[molecule].phi = 135.0*(3.141592653589/180.0);
                          coordinates[molecule].tetta = 70.0*(3.141592653589/180.0);
                          molecule++;
                         }
            }
            else
            {
           if((j%2)!=0){
                        coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                        coordinates[molecule].y = PBC2D(Ly, i*b);
                        coordinates[molecule].phi = -135.0*(3.141592653589/180.0);
                        coordinates[molecule].tetta = 70.0*(3.141592653589/180.0);
                        molecule++;
                       }
            }
        }
    }

 nPart = molecule;
 density = (166.113/3.318/3.318)*nPart/(Lx*Ly);

 cout << "Herringbone structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mmol/m2" << "\t" << "Lx/Ly: " << Lx/Ly << endl;
}
