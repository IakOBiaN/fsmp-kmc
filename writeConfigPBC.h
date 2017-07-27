  using namespace std;

  void writeConfigPBC (int &nPart, double &sigma, double &Lx, double &Ly, vector <state> &coordinates, double &write_rad, string config)
  {
   stringstream name;
   name << "xyz_d" << nPart/Lx/Ly << "_Lx" << Lx << "_Ly" << Ly << "_" << config << ".dat";
   ofstream fileOutput(name.str().c_str(), ios_base::app);

   double x_up, x_down, y_up, y_down, z_up, z_down;
   fileOutput << "x" << "\t" << "y" << "\t" << "z" << "\t" << "r" << endl;
   for(int i = 0; i < nPart; i++)
      {
        x_up = coordinates[i].x + (write_rad)*sin(coordinates[i].tetta)*cos(coordinates[i].phi);
        y_up = coordinates[i].y + (write_rad)*sin(coordinates[i].tetta)*sin(coordinates[i].phi);
        z_up = (write_rad)*cos(coordinates[i].tetta);
        fileOutput << x_up << "\t" << y_up << "\t" << z_up << "\t" << sigma/8 << endl;
        fileOutput << x_up << "\t" << y_up - Ly << "\t" << z_up << "\t" << sigma/8 << endl;
        fileOutput << x_up << "\t" << y_up + Ly << "\t" << z_up << "\t" << sigma/8 << endl;
        fileOutput << x_up - Lx << "\t" << y_up << "\t" << z_up << "\t" << sigma/8 << endl;
        fileOutput << x_up + Lx << "\t" << y_up << "\t" << z_up << "\t" << sigma/8 << endl;

        x_down = coordinates[i].x - (write_rad)*sin(coordinates[i].tetta)*cos(coordinates[i].phi);
        y_down = coordinates[i].y - (write_rad)*sin(coordinates[i].tetta)*sin(coordinates[i].phi);
        z_down = -(write_rad)*cos(coordinates[i].tetta);
        fileOutput << x_down << "\t" << y_down << "\t" << z_down << "\t" << sigma/8 << endl;
        fileOutput << x_down << "\t" << y_down - Ly << "\t" << z_down << "\t" << sigma/8 << endl;
        fileOutput << x_down << "\t" << y_down + Ly << "\t" << z_down << "\t" << sigma/8 << endl;
        fileOutput << x_down - Lx << "\t" << y_down << "\t" << z_down << "\t" << sigma/8 << endl;
        fileOutput << x_down + Lx << "\t" << y_down << "\t" << z_down << "\t" << sigma/8 << endl;
      }
   fileOutput.close();
  }
