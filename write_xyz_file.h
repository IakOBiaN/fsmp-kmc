  void write_xyz_file (int &nPart, double &Lx, double &Ly, double &temperature, vector <state> &coordinates, int frame, double distance, bool init)
{

   stringstream name;
   name << "xyz_N" << nPart << "_Lx=" << Lx << "_Ly=" << Ly << "_T=" << temperature << ".xyz";
   if (init) {ofstream fileOutput(name.str().c_str(), ios_base::trunc);fileOutput.close();}
   else {
        ofstream fileOutput(name.str().c_str(), ios_base::app);
        fileOutput << nPart*4 << endl;
        fileOutput << "Lattice=\"" << Lx << " 0 0 0 " << Ly << " 0 0 0 1\"" << endl;
        double dn2 = 3.536/2.0;

        for(int i = 0; i < nPart; i++)
          {
            double dop = coordinates[i].x + dn2 * coordinates[i].cos_phi;
            if (dop > Lx) {dop -= Lx;}
            if (dop < 0)  {dop += Lx;}
            double xxx = dop;
            dop = coordinates[i].y + dn2 * coordinates[i].sin_phi;
            if (dop > Ly) {dop -= Ly;}
            if (dop < 0)  {dop += Ly;}
            double yyy = dop;
            fileOutput << "N " << xxx << " " << yyy << " " << 0 << endl;

            dop = coordinates[i].x + dn2 * (-0.5*coordinates[i].cos_phi - 0.86602540378443864676372317075294*coordinates[i].sin_phi);
            if (dop > Lx) {dop -= Lx;}
            if (dop < 0)  {dop += Lx;}
            xxx = dop;
            dop = coordinates[i].y + dn2 * (-0.5*coordinates[i].sin_phi + 0.86602540378443864676372317075294*coordinates[i].cos_phi);
            if (dop > Ly) {dop -= Ly;}
            if (dop < 0)  {dop += Ly;}
            yyy = dop;
            fileOutput << "N " << xxx << " " << yyy << " " << 0 << endl;

            dop = coordinates[i].x + dn2 * (-0.5*coordinates[i].cos_phi + 0.86602540378443864676372317075294*coordinates[i].sin_phi);
            if (dop > Lx) {dop -= Lx;}
            if (dop < 0)  {dop += Lx;}
            xxx = dop;
            dop = coordinates[i].y + dn2 * (-0.5*coordinates[i].sin_phi - 0.86602540378443864676372317075294*coordinates[i].cos_phi);
            if (dop > Ly) {dop -= Ly;}
            if (dop < 0)  {dop += Ly;}
            yyy = dop;
            fileOutput << "N " << xxx << " " << yyy << " " << 0 << endl;

            fileOutput << "H " << coordinates[i].x << " " << coordinates[i].y << " " << 0 << endl;

          }
        fileOutput.close();
        }
}
