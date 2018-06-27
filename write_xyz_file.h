  void write_xyz_file_N2 (int &nPart, double &Lx, double &Ly, double &temperature, vector <state> &coordinates, int frame, double distance, bool init)
{

   stringstream name;
   name << "xyz_N" << nPart << "_Lx=" << Lx << "_Ly=" << Ly << "_T=" << temperature << ".xyz";
   if (init) {ofstream fileOutput(name.str().c_str(), ios_base::trunc);fileOutput.close();}
   else {
        ofstream fileOutput(name.str().c_str(), ios_base::app);
        fileOutput << nPart*2 << endl;
        fileOutput << frame << endl;
        fileOutput << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\"" << endl;

        for(int i = 0; i < nPart; i++)
          {
            fileOutput << "N " << coordinates[i].x+coordinates[i].cos_phi*distance/2 << " " << coordinates[i].y+coordinates[i].sin_phi*distance/2 << " " << 0 << endl;
            fileOutput << "N " << coordinates[i].x-coordinates[i].cos_phi*distance/2 << " " << coordinates[i].y-coordinates[i].sin_phi*distance/2 << " " << 0 << endl;
          }
        fileOutput.close();
        }
}
