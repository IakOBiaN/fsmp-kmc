using namespace std;

initConfigHerringbone (int &nPart, double &density, double &sigma, vector <state> &coordinates, double &Lx, double &Ly)
{

 // Get the corresponding box size
 //L = 2*nPart/density; // Lx*Ly
 //Ly = sqrt(L*sqrt(3)/2);
 //Lx = Ly*2/sqrt(3);

 // Generate the box size
 // based on the number of molecules
 //double LxLy = (15.72*nPart/3.318/3.318)*density;
 //Ly = sqrt(LxLy*sqrt(3)/2);
 //Lx = Ly*2/sqrt(3);

 //cout << "L: " << L << "\t" << "Lx: " << Lx << "\t" << "Ly: " << Ly << endl;

 double a = 4.04/3.318/2; // half of x-unit vector of herringbone unit cell in sigma
 double b = 7.00/3.318/2; // half of y-unit vector of herringbone unit cell in sigma
 double gamma = 88.5*(3.141592653589/180.0);
 //int stepsX = Lx/a;
 //int stepsY = Ly/b;
 int stepsX = 40;
 int stepsY = 20;
 Lx = a * 40;
 Ly = b * 20;

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
                         coordinates[molecule].phi = 135.0*(3.141592653589/180.0);
                         coordinates[molecule].tetta = 70.0*(3.141592653589/180.0);
                         molecule++;
                        }
           }
         else
          {
           if((j%2)==0){
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
 density = 15.72/(Lx*Ly*sin(gamma)*3.318*3.318/nPart);

 cout << "N: " << molecule << " density: " << density << endl;
}
