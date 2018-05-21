void write_xy_matrix(int &nPart, double &Lx, double &Ly, double &temperature, vector <vector <double> > &xy_matrix)
{
   stringstream name;
   name <<  "N" << nPart << "_Lx" << Lx << "_Ly" << Ly << "_T" << temperature << ".dat";
   ofstream fileOutput(name.str().c_str(), ios_base::trunc);

   int Xmax = Lx*5;
   int Ymax = Ly*5;

   for(int i = 0; i < Xmax; i++)
      {
          for (int j = 0; j < Ymax; j++)
          {
              fileOutput << xy_matrix[i][j] << " ";
          }
          fileOutput << endl;
      }
   fileOutput.close();
}
