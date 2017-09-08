 void virial_pressure (int &nPart, double &Lx, double &Ly, double &beta, double &Rc, double &Rc2, vector <state> &coordinates,
                         double &p_N, double &p_T, double &p_Tot, const double &A, double &C_q)
 {

  // Loop over all distinct particle pairs
  //double dr = 0.01;
  for (int molA = 0; molA < nPart-1; molA++)
  {
     for (int molB = molA+1; molB < nPart; molB++)
     {
         inter_for_pressure_ab(coordinates[molA], coordinates[molB], Rc, Rc2, Lx, Ly, A, C_q, beta, p_N, p_T, p_Tot);
         //double e_old = Inter_potential(coordinates[molA], coordinates[molB], Rc, Rc2, Lx, Ly, A, C_q, beta);
         //state trialPart = coordinates[molB];
         //double dist = sqrt(pow(coordinates[molA].x - coordinates[molB].x, 2) + pow(coordinates[molA].y - coordinates[molB].y, 2));
         //double dx = dr*(coordinates[molA].x - coordinates[molB].x)/dist;
         //trialPart.x -= dx;
         //double dy = dr*(coordinates[molA].y - coordinates[molB].y)/dist;
         //trialPart.y -= dy;
         //double e_new = Inter_potential(coordinates[molA], trialPart, Rc, Rc2, Lx, Ly, A, C_q, beta);
         //p_N += (e_new - e_old)*(coordinates[molA].x - coordinates[molB].x)/dist;
         //p_T += (e_new - e_old)*(coordinates[molA].y - coordinates[molB].y)/dist;
     }
  }
 //p_N = - p_N/dr/Lx/Ly;
 //p_T = - p_T/dr/Lx/Ly;
}
