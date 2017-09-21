using namespace std;

void initConfigPinwheel (int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double &coeff)
{

 double a = coeff*8.18/3.318/2; // half of x-unit vector of herringbone unit cell in sigma // ideal = 0.933
 double b = a*sqrt(3)/2.0; // half of y-unit vector of herringbone unit cell in sigma
 double gamma = 60.0*(3.141592653589/180.0);

 int nUnitCells = sqrt(nPart/4.0);

 Lx = nUnitCells*a*2.0;
 Ly = nUnitCells*b*2.0;

 int stepsX = floor(Lx/a + 0.5);
 int stepsY = floor(Ly/b + 0.5);

 cout << "Lx: " << Lx << "\t" << "Ly: " << Ly << endl;

 int molecule = 0; // Molecules counter

 for(int i = 0; i < stepsY; i++)
    {
     for(int j = 0; j < stepsX; j++)
        {
         if((i%2)==0)
           {
            if((j%2)!=0){
                         coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                         coordinates[molecule].y = PBC2D(Ly, i*b);
                         coordinates[molecule].phi = -150.0*(3.141592653589/180.0);
                         molecule++;
                        }
                    else{
                         coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                         coordinates[molecule].y = PBC2D(Ly, i*b);
                         coordinates[molecule].phi = 0.0;
                         molecule++;
                        }
           }
         else
          {
           if((j%2)!=0){
                        coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                        coordinates[molecule].y = PBC2D(Ly, i*b);
                        coordinates[molecule].phi = 90.0*(3.141592653589/180.0);
                        molecule++;
                       }
                   else{
                         coordinates[molecule].x = PBC2D(Lx, j*a - i*b*cos(gamma));
                         coordinates[molecule].y = PBC2D(Ly, i*b);
                         coordinates[molecule].phi = -30.0*(3.141592653589/180.0);
                         molecule++;
                        }
          }
        }
    }
 nPart = molecule;
 density = (166.113/3.318/3.318)*nPart/(Lx*Ly);

 cout << "Pinwheel structure: " << endl;
 cout << "N: " << molecule << "\t" << "density: " << density << " mmol/m2" << "\t" << "Lx/Ly: " << Lx/Ly << endl;
}
