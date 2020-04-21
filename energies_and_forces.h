double energy_calculation(state &molA, state &molB, double &r, double &ang1, double &ang2)
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
  double d_charges = 4.2; // A
  double QQ_const = Ke*q2*N_a;

  /// Lennard-Johns interaction ///
  double r_r0 = r - r0;
  double sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
  double r_r0_6 = r_r0*r_r0*r_r0*r_r0*r_r0*r_r0;
  double lj_dist_6 = sigma_r0_6/r_r0_6;
  U_LJ += 4*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol

  /// Coulomb's interaction ///
  double sin_half_carbox = sin(alpha/180.0*PI);
  double cos_half_carbox = cos(alpha/180.0*PI);

  // Coordinates of charges in the molecule A
  double sin_ang1 = sin(ang1/180.0*PI);
  double cos_ang1 = cos(ang1/180.0*PI);
  // for positive charges
  double molA_sin_add_half_carbox_p = sin_ang1*cos_half_carbox + cos_ang1*sin_half_carbox;
  double molA_cos_add_half_carbox_p = cos_ang1*cos_half_carbox - sin_ang1*sin_half_carbox;
  double molA_sin_add_half_carbox_p_120 = molA_sin_add_half_carbox_p*cos(120/180.0*PI) + molA_cos_add_half_carbox_p*sin(120/180.0*PI);
  double molA_cos_add_half_carbox_p_120 = molA_cos_add_half_carbox_p*cos(120/180.0*PI) - molA_sin_add_half_carbox_p*sin(120/180.0*PI);
  double molA_sin_add_half_carbox_p_240 = molA_sin_add_half_carbox_p*cos(240/180.0*PI) + molA_cos_add_half_carbox_p*sin(240/180.0*PI);
  double molA_cos_add_half_carbox_p_240 = molA_cos_add_half_carbox_p*cos(240/180.0*PI) - molA_sin_add_half_carbox_p*sin(240/180.0*PI);

  // for negative charges
  double molA_sin_add_half_carbox_n = sin_ang1*cos_half_carbox - cos_ang1*sin_half_carbox;
  double molA_cos_add_half_carbox_n = cos_ang1*cos_half_carbox + sin_ang1*sin_half_carbox;
  double molA_sin_add_half_carbox_n_120 = molA_sin_add_half_carbox_n*cos(120/180.0*PI) + molA_cos_add_half_carbox_n*sin(120/180.0*PI);
  double molA_cos_add_half_carbox_n_120 = molA_cos_add_half_carbox_n*cos(120/180.0*PI) - molA_sin_add_half_carbox_n*sin(120/180.0*PI);
  double molA_sin_add_half_carbox_n_240 = molA_sin_add_half_carbox_n*cos(240/180.0*PI) + molA_cos_add_half_carbox_n*sin(240/180.0*PI);
  double molA_cos_add_half_carbox_n_240 = molA_cos_add_half_carbox_n*cos(240/180.0*PI) - molA_sin_add_half_carbox_n*sin(240/180.0*PI);

  // x,y - coordinates of six charges (3 positive and 3 negative) in the A molecule
  double molA_x_q1_p = molA.x + d_charges*molA_cos_add_half_carbox_p;
  double molA_y_q1_p = molA.y + d_charges*molA_sin_add_half_carbox_p;
  double molA_x_q2_p = molA.x + d_charges*molA_cos_add_half_carbox_p_120;
  double molA_y_q2_p = molA.y + d_charges*molA_sin_add_half_carbox_p_120;
  double molA_x_q3_p = molA.x + d_charges*molA_cos_add_half_carbox_p_240;
  double molA_y_q3_p = molA.y + d_charges*molA_sin_add_half_carbox_p_240;
  double molA_x_q1_n = molA.x + d_charges*molA_cos_add_half_carbox_n;
  double molA_y_q1_n = molA.y + d_charges*molA_sin_add_half_carbox_n;
  double molA_x_q2_n = molA.x + d_charges*molA_cos_add_half_carbox_n_120;
  double molA_y_q2_n = molA.y + d_charges*molA_sin_add_half_carbox_n_120;
  double molA_x_q3_n = molA.x + d_charges*molA_cos_add_half_carbox_n_240;
  double molA_y_q3_n = molA.y + d_charges*molA_sin_add_half_carbox_n_240;

  // Coordinates of charges in the molecule B
  double sin_ang2 = sin(ang2/180.0*PI);
  double cos_ang2 = cos(ang2/180.0*PI);
  // for positive charges
  double molB_sin_add_half_carbox_p = sin_ang2*cos_half_carbox + cos_ang2*sin_half_carbox;
  double molB_cos_add_half_carbox_p = cos_ang2*cos_half_carbox - sin_ang2*sin_half_carbox;
  double molB_sin_add_half_carbox_p_120 = molB_sin_add_half_carbox_p*cos(120/180.0*PI) + molB_cos_add_half_carbox_p*sin(120/180.0*PI);
  double molB_cos_add_half_carbox_p_120 = molB_cos_add_half_carbox_p*cos(120/180.0*PI) - molB_sin_add_half_carbox_p*sin(120/180.0*PI);
  double molB_sin_add_half_carbox_p_240 = molB_sin_add_half_carbox_p*cos(240/180.0*PI) + molB_cos_add_half_carbox_p*sin(240/180.0*PI);
  double molB_cos_add_half_carbox_p_240 = molB_cos_add_half_carbox_p*cos(240/180.0*PI) - molB_sin_add_half_carbox_p*sin(240/180.0*PI);
  // for negative charges
  double molB_sin_add_half_carbox_n = sin_ang2*cos_half_carbox - cos_ang2*sin_half_carbox;
  double molB_cos_add_half_carbox_n = cos_ang2*cos_half_carbox + sin_ang2*sin_half_carbox;
  double molB_sin_add_half_carbox_n_120 = molB_sin_add_half_carbox_n*cos(120/180.0*PI) + molB_cos_add_half_carbox_n*sin(120/180.0*PI);
  double molB_cos_add_half_carbox_n_120 = molB_cos_add_half_carbox_n*cos(120/180.0*PI) - molB_sin_add_half_carbox_n*sin(120/180.0*PI);
  double molB_sin_add_half_carbox_n_240 = molB_sin_add_half_carbox_n*cos(240/180.0*PI) + molB_cos_add_half_carbox_n*sin(240/180.0*PI);
  double molB_cos_add_half_carbox_n_240 = molB_cos_add_half_carbox_n*cos(240/180.0*PI) - molB_sin_add_half_carbox_n*sin(240/180.0*PI);
  // x,y - coordinates of six charges (3 positive and 3 negative) in the B molecule
  double molB_x_q1_p = molB.x + d_charges*molB_cos_add_half_carbox_p;
  double molB_y_q1_p = molB.y + d_charges*molB_sin_add_half_carbox_p;
  double molB_x_q2_p = molB.x + d_charges*molB_cos_add_half_carbox_p_120;
  double molB_y_q2_p = molB.y + d_charges*molB_sin_add_half_carbox_p_120;
  double molB_x_q3_p = molB.x + d_charges*molB_cos_add_half_carbox_p_240;
  double molB_y_q3_p = molB.y + d_charges*molB_sin_add_half_carbox_p_240;
  double molB_x_q1_n = molB.x + d_charges*molB_cos_add_half_carbox_n;
  double molB_y_q1_n = molB.y + d_charges*molB_sin_add_half_carbox_n;
  double molB_x_q2_n = molB.x + d_charges*molB_cos_add_half_carbox_n_120;
  double molB_y_q2_n = molB.y + d_charges*molB_sin_add_half_carbox_n_120;
  double molB_x_q3_n = molB.x + d_charges*molB_cos_add_half_carbox_n_240;
  double molB_y_q3_n = molB.y + d_charges*molB_sin_add_half_carbox_n_240;

  // Interaction energy between 6 charges of the A molecule and 6 charges of the B molecules
  // Positive charges of A molecule and positive charges of B molecule
  U_QQ += QQ_const/sqrt((molB_x_q1_p - molA_x_q1_p)*(molB_x_q1_p - molA_x_q1_p) + (molB_y_q1_p - molA_y_q1_p)*(molB_y_q1_p - molA_y_q1_p));
  U_QQ += QQ_const/sqrt((molB_x_q2_p - molA_x_q1_p)*(molB_x_q2_p - molA_x_q1_p) + (molB_y_q2_p - molA_y_q1_p)*(molB_y_q2_p - molA_y_q1_p));
  U_QQ += QQ_const/sqrt((molB_x_q3_p - molA_x_q1_p)*(molB_x_q3_p - molA_x_q1_p) + (molB_y_q3_p - molA_y_q1_p)*(molB_y_q3_p - molA_y_q1_p));
  U_QQ += QQ_const/sqrt((molB_x_q1_p - molA_x_q2_p)*(molB_x_q1_p - molA_x_q2_p) + (molB_y_q1_p - molA_y_q2_p)*(molB_y_q1_p - molA_y_q2_p));
  U_QQ += QQ_const/sqrt((molB_x_q2_p - molA_x_q2_p)*(molB_x_q2_p - molA_x_q2_p) + (molB_y_q2_p - molA_y_q2_p)*(molB_y_q2_p - molA_y_q2_p));
  U_QQ += QQ_const/sqrt((molB_x_q3_p - molA_x_q2_p)*(molB_x_q3_p - molA_x_q2_p) + (molB_y_q3_p - molA_y_q2_p)*(molB_y_q3_p - molA_y_q2_p));
  U_QQ += QQ_const/sqrt((molB_x_q1_p - molA_x_q3_p)*(molB_x_q1_p - molA_x_q3_p) + (molB_y_q1_p - molA_y_q3_p)*(molB_y_q1_p - molA_y_q3_p));
  U_QQ += QQ_const/sqrt((molB_x_q2_p - molA_x_q3_p)*(molB_x_q2_p - molA_x_q3_p) + (molB_y_q2_p - molA_y_q3_p)*(molB_y_q2_p - molA_y_q3_p));
  U_QQ += QQ_const/sqrt((molB_x_q3_p - molA_x_q3_p)*(molB_x_q3_p - molA_x_q3_p) + (molB_y_q3_p - molA_y_q3_p)*(molB_y_q3_p - molA_y_q3_p));
  // Negative charges of A molecule and negative charges of B molecule
  U_QQ += QQ_const/sqrt((molB_x_q1_n - molA_x_q1_n)*(molB_x_q1_n - molA_x_q1_n) + (molB_y_q1_n - molA_y_q1_n)*(molB_y_q1_n - molA_y_q1_n));
  U_QQ += QQ_const/sqrt((molB_x_q2_n - molA_x_q1_n)*(molB_x_q2_n - molA_x_q1_n) + (molB_y_q2_n - molA_y_q1_n)*(molB_y_q2_n - molA_y_q1_n));
  U_QQ += QQ_const/sqrt((molB_x_q3_n - molA_x_q1_n)*(molB_x_q3_n - molA_x_q1_n) + (molB_y_q3_n - molA_y_q1_n)*(molB_y_q3_n - molA_y_q1_n));
  U_QQ += QQ_const/sqrt((molB_x_q1_n - molA_x_q2_n)*(molB_x_q1_n - molA_x_q2_n) + (molB_y_q1_n - molA_y_q2_n)*(molB_y_q1_n - molA_y_q2_n));
  U_QQ += QQ_const/sqrt((molB_x_q2_n - molA_x_q2_n)*(molB_x_q2_n - molA_x_q2_n) + (molB_y_q2_n - molA_y_q2_n)*(molB_y_q2_n - molA_y_q2_n));
  U_QQ += QQ_const/sqrt((molB_x_q3_n - molA_x_q2_n)*(molB_x_q3_n - molA_x_q2_n) + (molB_y_q3_n - molA_y_q2_n)*(molB_y_q3_n - molA_y_q2_n));
  U_QQ += QQ_const/sqrt((molB_x_q1_n - molA_x_q3_n)*(molB_x_q1_n - molA_x_q3_n) + (molB_y_q1_n - molA_y_q3_n)*(molB_y_q1_n - molA_y_q3_n));
  U_QQ += QQ_const/sqrt((molB_x_q2_n - molA_x_q3_n)*(molB_x_q2_n - molA_x_q3_n) + (molB_y_q2_n - molA_y_q3_n)*(molB_y_q2_n - molA_y_q3_n));
  U_QQ += QQ_const/sqrt((molB_x_q3_n - molA_x_q3_n)*(molB_x_q3_n - molA_x_q3_n) + (molB_y_q3_n - molA_y_q3_n)*(molB_y_q3_n - molA_y_q3_n));
  // Negative charges of A molecule and positive charges of B molecule
  U_QQ += -QQ_const/sqrt((molB_x_q1_p - molA_x_q1_n)*(molB_x_q1_p - molA_x_q1_n) + (molB_y_q1_p - molA_y_q1_n)*(molB_y_q1_p - molA_y_q1_n));
  U_QQ += -QQ_const/sqrt((molB_x_q2_p - molA_x_q1_n)*(molB_x_q2_p - molA_x_q1_n) + (molB_y_q2_p - molA_y_q1_n)*(molB_y_q2_p - molA_y_q1_n));
  U_QQ += -QQ_const/sqrt((molB_x_q3_p - molA_x_q1_n)*(molB_x_q3_p - molA_x_q1_n) + (molB_y_q3_p - molA_y_q1_n)*(molB_y_q3_p - molA_y_q1_n));
  U_QQ += -QQ_const/sqrt((molB_x_q1_p - molA_x_q2_n)*(molB_x_q1_p - molA_x_q2_n) + (molB_y_q1_p - molA_y_q2_n)*(molB_y_q1_p - molA_y_q2_n));
  U_QQ += -QQ_const/sqrt((molB_x_q2_p - molA_x_q2_n)*(molB_x_q2_p - molA_x_q2_n) + (molB_y_q2_p - molA_y_q2_n)*(molB_y_q2_p - molA_y_q2_n));
  U_QQ += -QQ_const/sqrt((molB_x_q3_p - molA_x_q2_n)*(molB_x_q3_p - molA_x_q2_n) + (molB_y_q3_p - molA_y_q2_n)*(molB_y_q3_p - molA_y_q2_n));
  U_QQ += -QQ_const/sqrt((molB_x_q1_p - molA_x_q3_n)*(molB_x_q1_p - molA_x_q3_n) + (molB_y_q1_p - molA_y_q3_n)*(molB_y_q1_p - molA_y_q3_n));
  U_QQ += -QQ_const/sqrt((molB_x_q2_p - molA_x_q3_n)*(molB_x_q2_p - molA_x_q3_n) + (molB_y_q2_p - molA_y_q3_n)*(molB_y_q2_p - molA_y_q3_n));
  U_QQ += -QQ_const/sqrt((molB_x_q3_p - molA_x_q3_n)*(molB_x_q3_p - molA_x_q3_n) + (molB_y_q3_p - molA_y_q3_n)*(molB_y_q3_p - molA_y_q3_n));
  // Positive charges of A molecule and negative charges of B molecule
  U_QQ += -QQ_const/sqrt((molB_x_q1_n - molA_x_q1_p)*(molB_x_q1_n - molA_x_q1_p) + (molB_y_q1_n - molA_y_q1_p)*(molB_y_q1_n - molA_y_q1_p));
  U_QQ += -QQ_const/sqrt((molB_x_q2_n - molA_x_q1_p)*(molB_x_q2_n - molA_x_q1_p) + (molB_y_q2_n - molA_y_q1_p)*(molB_y_q2_n - molA_y_q1_p));
  U_QQ += -QQ_const/sqrt((molB_x_q3_n - molA_x_q1_p)*(molB_x_q3_n - molA_x_q1_p) + (molB_y_q3_n - molA_y_q1_p)*(molB_y_q3_n - molA_y_q1_p));
  U_QQ += -QQ_const/sqrt((molB_x_q1_n - molA_x_q2_p)*(molB_x_q1_n - molA_x_q2_p) + (molB_y_q1_n - molA_y_q2_p)*(molB_y_q1_n - molA_y_q2_p));
  U_QQ += -QQ_const/sqrt((molB_x_q2_n - molA_x_q2_p)*(molB_x_q2_n - molA_x_q2_p) + (molB_y_q2_n - molA_y_q2_p)*(molB_y_q2_n - molA_y_q2_p));
  U_QQ += -QQ_const/sqrt((molB_x_q3_n - molA_x_q2_p)*(molB_x_q3_n - molA_x_q2_p) + (molB_y_q3_n - molA_y_q2_p)*(molB_y_q3_n - molA_y_q2_p));
  U_QQ += -QQ_const/sqrt((molB_x_q1_n - molA_x_q3_p)*(molB_x_q1_n - molA_x_q3_p) + (molB_y_q1_n - molA_y_q3_p)*(molB_y_q1_n - molA_y_q3_p));
  U_QQ += -QQ_const/sqrt((molB_x_q2_n - molA_x_q3_p)*(molB_x_q2_n - molA_x_q3_p) + (molB_y_q2_n - molA_y_q3_p)*(molB_y_q2_n - molA_y_q3_p));
  U_QQ += -QQ_const/sqrt((molB_x_q3_n - molA_x_q3_p)*(molB_x_q3_n - molA_x_q3_p) + (molB_y_q3_n - molA_y_q3_p)*(molB_y_q3_n - molA_y_q3_p));

  return U_LJ + U_QQ;
}


results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta, bool pressure_calculation)
{

    double max_dist = 18.0;
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
                 double dang = dist_x/r;	// Calculate the cosine of the angle between OX and distance vector
                 if (dist_y<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
                 ang1 = ang_molA - dang; // Angle of A molecule when distance vector is parallel to the OX axis
                 ang2 = ang_molB - dang; // Angle of A molecule when distance vector is parallel to the OX axis
                 // Molecule should always has the angle in the range of 0-360 degrees
                 if (ang1<0) {ang1 += 360.0;}
                 if (ang2<0) {ang2 += 360.0;}
                 if (ang1>359.999) {ang1 -= 360.0;}
                 if (ang2>359.999) {ang2 -= 360.0;}

                 energy = energy + energy_calculation(molA, molB, r, ang1, ang2);
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
                       double dang = dist_x_plus_delta/r;
                       if (dist_y<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
                       ang1 = ang_molA - dang;
                       ang2 = ang_molB - dang;
                       if (ang1<0) {ang1 += 360.0;}
                       if (ang2<0) {ang2 += 360.0;}
                       if (ang1>359.5) {ang1 -= 360.0;}
                       if (ang2>359.5) {ang2 -= 360.0;}
                       var_press[numerator] = energy_calculation(molA, molB, r, ang1, ang2);
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
                        double dang = dist_x/r;
                        if (dist_y_plus_delta<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
                        ang1 = ang_molA - dang;
                        ang2 = ang_molB - dang;
                        if (ang1<0) {ang1 += 360.0;}
                        if (ang2<0) {ang2 += 360.0;}
                        if (ang1>359.5) {ang1 -= 360.0;}
                        if (ang2>359.5) {ang2 -= 360.0;}
                        var_press[numerator] = energy_calculation(molA, molB, r, ang1, ang2);
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
