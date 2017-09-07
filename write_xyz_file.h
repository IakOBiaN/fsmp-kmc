  void write_xyz_file (int &nPart, double &Lx, double &Ly, double &temperature, vector <state> &coordinates, int frame, double distance, bool init)
{

   stringstream name;
   name << "xyz_N" << nPart << "_Lx=" << Lx << "_Ly=" << Ly << "_T=" << temperature << ".xyz";
   if (init) {ofstream fileOutput(name.str().c_str(), ios_base::trunc);fileOutput.close();}
   else {
        ofstream fileOutput(name.str().c_str(), ios_base::app);
        fileOutput << nPart*2 << endl;
        fileOutput << frame << endl;
        for(int i = 0; i < nPart; i++)
          {
            if (abs(sin(coordinates[i].tetta))<0.7)
            {
            fileOutput << "O " << coordinates[i].x+sin(coordinates[i].tetta)*cos(coordinates[i].phi)*distance/2 << " " << coordinates[i].y+sin(coordinates[i].tetta)*sin(coordinates[i].phi)*distance/2 << " " << cos(coordinates[i].tetta)*distance/2 << endl;
            fileOutput << "O " << coordinates[i].x-sin(coordinates[i].tetta)*cos(coordinates[i].phi)*distance/2 << " " << coordinates[i].y-sin(coordinates[i].tetta)*sin(coordinates[i].phi)*distance/2 << " " << -cos(coordinates[i].tetta)*distance/2 << endl;
            }
            else
            {
            fileOutput << "N " << coordinates[i].x+sin(coordinates[i].tetta)*cos(coordinates[i].phi)*distance/2 << " " << coordinates[i].y+sin(coordinates[i].tetta)*sin(coordinates[i].phi)*distance/2 << " " << cos(coordinates[i].tetta)*distance/2 << endl;
            fileOutput << "N " << coordinates[i].x-sin(coordinates[i].tetta)*cos(coordinates[i].phi)*distance/2 << " " << coordinates[i].y-sin(coordinates[i].tetta)*sin(coordinates[i].phi)*distance/2 << " " << -cos(coordinates[i].tetta)*distance/2 << endl;
            }
          }
       fileOutput.close();
       }
}
