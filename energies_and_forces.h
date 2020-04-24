double energy_calculation(state &molA, state &molB, double &r)
{
  double U_LJ = 0, U_QQ = 0;
  double alpha = 15.0; // Half of angle between oxygens in the carboxylic group of TMA molecule
  double gamma = 120.0;
  double double_gamma = 240.0;
  double r0 = 7.5887; // Hard core radius in A
  double sigma_r0 = 11.052 - 7.5887; // in A
  double eps = 82.525; // in J/mol or eps/k_B = 360 K
  const double Ke = 8.987551e+19; // in J*A/C^2
  double N_a = 6.02214e+23;
  double q = 0.45*1.6e-19; // in C
  double q2 = q*q;
  double QQ_const = Ke*q2*N_a;

  /// Lennard-Johns interaction ///
  double r_r0 = r - r0;
  double sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
  double r_r0_6 = r_r0*r_r0*r_r0*r_r0*r_r0*r_r0;
  double lj_dist_6 = sigma_r0_6/r_r0_6;
  U_LJ += 4*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol

  /// Coulomb's interaction ///

  // Interaction energy between 6 charges of the A molecule and 6 charges of the B molecules
  // Positive charges of A molecule and positive charges of B molecule
  U_QQ += QQ_const/sqrt((molB.q1x_p - molA.q1x_p)*(molB.q1x_p - molA.q1x_p) + (molB.q1y_p - molA.q1y_p)*(molB.q1y_p - molA.q1y_p));
  U_QQ += QQ_const/sqrt((molB.q2x_p - molA.q1x_p)*(molB.q2x_p - molA.q1x_p) + (molB.q2y_p - molA.q1y_p)*(molB.q2y_p - molA.q1y_p));
  U_QQ += QQ_const/sqrt((molB.q3x_p - molA.q1x_p)*(molB.q3x_p - molA.q1x_p) + (molB.q3y_p - molA.q1y_p)*(molB.q3y_p - molA.q1y_p));
  U_QQ += QQ_const/sqrt((molB.q1x_p - molA.q2x_p)*(molB.q1x_p - molA.q2x_p) + (molB.q1y_p - molA.q2y_p)*(molB.q1y_p - molA.q2y_p));
  U_QQ += QQ_const/sqrt((molB.q2x_p - molA.q2x_p)*(molB.q2x_p - molA.q2x_p) + (molB.q2y_p - molA.q2y_p)*(molB.q2y_p - molA.q2y_p));
  U_QQ += QQ_const/sqrt((molB.q3x_p - molA.q2x_p)*(molB.q3x_p - molA.q2x_p) + (molB.q3y_p - molA.q2y_p)*(molB.q3y_p - molA.q2y_p));
  U_QQ += QQ_const/sqrt((molB.q1x_p - molA.q3x_p)*(molB.q1x_p - molA.q3x_p) + (molB.q1y_p - molA.q3y_p)*(molB.q1y_p - molA.q3y_p));
  U_QQ += QQ_const/sqrt((molB.q2x_p - molA.q3x_p)*(molB.q2x_p - molA.q3x_p) + (molB.q2y_p - molA.q3y_p)*(molB.q2y_p - molA.q3y_p));
  U_QQ += QQ_const/sqrt((molB.q3x_p - molA.q3x_p)*(molB.q3x_p - molA.q3x_p) + (molB.q3y_p - molA.q3y_p)*(molB.q3y_p - molA.q3y_p));
  // Negative charges of A molecule and negative charges of B molecule
  U_QQ += QQ_const/sqrt((molB.q1x_n - molA.q1x_n)*(molB.q1x_n - molA.q1x_n) + (molB.q1y_n - molA.q1y_n)*(molB.q1y_n - molA.q1y_n));
  U_QQ += QQ_const/sqrt((molB.q2x_n - molA.q1x_n)*(molB.q2x_n - molA.q1x_n) + (molB.q2y_n - molA.q1y_n)*(molB.q2y_n - molA.q1y_n));
  U_QQ += QQ_const/sqrt((molB.q3x_n - molA.q1x_n)*(molB.q3x_n - molA.q1x_n) + (molB.q3y_n - molA.q1y_n)*(molB.q3y_n - molA.q1y_n));
  U_QQ += QQ_const/sqrt((molB.q1x_n - molA.q2x_n)*(molB.q1x_n - molA.q2x_n) + (molB.q1y_n - molA.q2y_n)*(molB.q1y_n - molA.q2y_n));
  U_QQ += QQ_const/sqrt((molB.q2x_n - molA.q2x_n)*(molB.q2x_n - molA.q2x_n) + (molB.q2y_n - molA.q2y_n)*(molB.q2y_n - molA.q2y_n));
  U_QQ += QQ_const/sqrt((molB.q3x_n - molA.q2x_n)*(molB.q3x_n - molA.q2x_n) + (molB.q3y_n - molA.q2y_n)*(molB.q3y_n - molA.q2y_n));
  U_QQ += QQ_const/sqrt((molB.q1x_n - molA.q3x_n)*(molB.q1x_n - molA.q3x_n) + (molB.q1y_n - molA.q3y_n)*(molB.q1y_n - molA.q3y_n));
  U_QQ += QQ_const/sqrt((molB.q2x_n - molA.q3x_n)*(molB.q2x_n - molA.q3x_n) + (molB.q2y_n - molA.q3y_n)*(molB.q2y_n - molA.q3y_n));
  U_QQ += QQ_const/sqrt((molB.q3x_n - molA.q3x_n)*(molB.q3x_n - molA.q3x_n) + (molB.q3y_n - molA.q3y_n)*(molB.q3y_n - molA.q3y_n));
  // Negative charges of A molecule and positive charges of B molecule
  U_QQ += -QQ_const/sqrt((molB.q1x_p - molA.q1x_n)*(molB.q1x_p - molA.q1x_n) + (molB.q1y_p - molA.q1y_n)*(molB.q1y_p - molA.q1y_n));
  U_QQ += -QQ_const/sqrt((molB.q2x_p - molA.q1x_n)*(molB.q2x_p - molA.q1x_n) + (molB.q2y_p - molA.q1y_n)*(molB.q2y_p - molA.q1y_n));
  U_QQ += -QQ_const/sqrt((molB.q3x_p - molA.q1x_n)*(molB.q3x_p - molA.q1x_n) + (molB.q3y_p - molA.q1y_n)*(molB.q3y_p - molA.q1y_n));
  U_QQ += -QQ_const/sqrt((molB.q1x_p - molA.q2x_n)*(molB.q1x_p - molA.q2x_n) + (molB.q1y_p - molA.q2y_n)*(molB.q1y_p - molA.q2y_n));
  U_QQ += -QQ_const/sqrt((molB.q2x_p - molA.q2x_n)*(molB.q2x_p - molA.q2x_n) + (molB.q2y_p - molA.q2y_n)*(molB.q2y_p - molA.q2y_n));
  U_QQ += -QQ_const/sqrt((molB.q3x_p - molA.q2x_n)*(molB.q3x_p - molA.q2x_n) + (molB.q3y_p - molA.q2y_n)*(molB.q3y_p - molA.q2y_n));
  U_QQ += -QQ_const/sqrt((molB.q1x_p - molA.q3x_n)*(molB.q1x_p - molA.q3x_n) + (molB.q1y_p - molA.q3y_n)*(molB.q1y_p - molA.q3y_n));
  U_QQ += -QQ_const/sqrt((molB.q2x_p - molA.q3x_n)*(molB.q2x_p - molA.q3x_n) + (molB.q2y_p - molA.q3y_n)*(molB.q2y_p - molA.q3y_n));
  U_QQ += -QQ_const/sqrt((molB.q3x_p - molA.q3x_n)*(molB.q3x_p - molA.q3x_n) + (molB.q3y_p - molA.q3y_n)*(molB.q3y_p - molA.q3y_n));
  // Positive charges of A molecule and negative charges of B molecule
  U_QQ += -QQ_const/sqrt((molB.q1x_n - molA.q1x_p)*(molB.q1x_n - molA.q1x_p) + (molB.q1y_n - molA.q1y_p)*(molB.q1y_n - molA.q1y_p));
  U_QQ += -QQ_const/sqrt((molB.q2x_n - molA.q1x_p)*(molB.q2x_n - molA.q1x_p) + (molB.q2y_n - molA.q1y_p)*(molB.q2y_n - molA.q1y_p));
  U_QQ += -QQ_const/sqrt((molB.q3x_n - molA.q1x_p)*(molB.q3x_n - molA.q1x_p) + (molB.q3y_n - molA.q1y_p)*(molB.q3y_n - molA.q1y_p));
  U_QQ += -QQ_const/sqrt((molB.q1x_n - molA.q2x_p)*(molB.q1x_n - molA.q2x_p) + (molB.q1y_n - molA.q2y_p)*(molB.q1y_n - molA.q2y_p));
  U_QQ += -QQ_const/sqrt((molB.q2x_n - molA.q2x_p)*(molB.q2x_n - molA.q2x_p) + (molB.q2y_n - molA.q2y_p)*(molB.q2y_n - molA.q2y_p));
  U_QQ += -QQ_const/sqrt((molB.q3x_n - molA.q2x_p)*(molB.q3x_n - molA.q2x_p) + (molB.q3y_n - molA.q2y_p)*(molB.q3y_n - molA.q2y_p));
  U_QQ += -QQ_const/sqrt((molB.q1x_n - molA.q3x_p)*(molB.q1x_n - molA.q3x_p) + (molB.q1y_n - molA.q3y_p)*(molB.q1y_n - molA.q3y_p));
  U_QQ += -QQ_const/sqrt((molB.q2x_n - molA.q3x_p)*(molB.q2x_n - molA.q3x_p) + (molB.q2y_n - molA.q3y_p)*(molB.q2y_n - molA.q3y_p));
  U_QQ += -QQ_const/sqrt((molB.q3x_n - molA.q3x_p)*(molB.q3x_n - molA.q3x_p) + (molB.q3y_n - molA.q3y_p)*(molB.q3y_n - molA.q3y_p));

  return U_LJ + U_QQ;
}


results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta, bool pressure_calculation)
{
    results en_and_press;
  	double ang1,ang2;
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
  	double energy = 0,pressure_X = 0, pressure_Y = 0;
  	double dist_x_plus_delta, dist_y_plus_delta; // Distance between A and B molecules including derivation step
    double delta = 0.1;
    double ang_molA = molA.phi;
    double ang_molB = molB.phi;
    double var_press[2];
    bool press_calc;
    int numerator;
    double U_all = 0;

    for (int id = -1; id < 2; id++)
    {
       dist_x = (molB.x - molA.x + id*Lx);
       if (abs(dist_x) > max_dist) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          dist_y = (molB.y - molA.y + jd*Ly);
          if (abs(dist_y) > max_dist) {continue;}
             r2 = dist_x*dist_x + dist_y*dist_y;
             r = sqrt(r2);
             if (r <= max_dist)
             {
                 energy = energy + energy_calculation(molA, molB, r);
             }
             if (pressure_calculation)
             {
                 numerator = 0;
                 press_calc = true;
                 for (int diff_x = -1; diff_x < 3; diff_x += 2)
                 {
                       dist_x_plus_delta = dist_x + diff_x*delta*2.0;
                       r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
                       r = sqrt(r2);
                       if (r > max_dist || r < min_dist) {press_calc=false;break;}

                       var_press[numerator] = energy_calculation(molA, molB, r);
                       numerator++;
                  }
                  if (press_calc) {pressure_X += (var_press[0]-var_press[1])/(delta*4.0)*dist_x;}
                  numerator = 0;
                  press_calc = true;
                  for (int diff_y = -1; diff_y < 3; diff_y += 2)
                  {
                        dist_y_plus_delta = dist_y + diff_y*delta*2.0;
                        r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
                        r = sqrt(r2);
                        if (r > max_dist || r < min_dist) {press_calc=false;break;}

                        var_press[numerator] = energy_calculation(molA, molB, r);
                        numerator++;
                   }
                   if (press_calc) {pressure_Y += (var_press[0]-var_press[1])/(delta*4.0)*dist_y;}
             }
       }
    }
    en_and_press.energy = energy;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
