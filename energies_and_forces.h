results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
    double ang1,ang2;
    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
		double r;	// Distance between A and B molecules
		double r2;	// Distance sqaured
    double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
		double energy = 0, var_energy = 0, pressure_X = 0, pressure_Y = 0;
    double dist_x_plus_delta, dist_y_plus_delta; // Distance between A and B molecules including derivation step
		double delta = dr/2.0;	// Derivation step is associated with distance step in numerical potential
		double var_press[2]; // Current energy values in numerical differenciation procedure
    int numerator;
    bool press_calc; // The pressure should  been calculated only when current distances lay in the working interval of num. potential
    double ang_molA = molA.phi;
    double ang_molB = molB.phi;
    double dist_n; // Float index in the numerical potential

    int dist,a1,a2; // indexes in the numerical potential array


    for (int id = -1; id < 2; id++)
    {
       dist_x = (molB_x - molA_x + id*Lx);
       if (abs(dist_x) > max_dist) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          dist_y = (molB_y - molA_y + jd*Ly);
          if (abs(dist_y) > max_dist) {continue;}
             r2 = dist_x*dist_x + dist_y*dist_y;
             r = sqrt(r2);
             if (r <= max_dist)
             {
               //calculation of energy for two molecules
               double dang = dist_x/r;	// Calculate the cosine of the angle between OX and distance vector
               if (dist_y<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
               ang1 = ang_molA - dang;
               ang2 = ang_molB - dang;
							 // Molecule should always has the angle in the range of 0-360 degrees
							 if (ang1<0) {ang1 += 360.0;}
               if (ang2<0) {ang2 += 360.0;}
               if (ang1>359.5) {ang1 -= 360.0;}
               if (ang2>359.5) {ang2 -= 360.0;}
               dist_n = (r-min_dist)/dr;
               dist = (int)(dist_n+0.5);
               a1 = (int)((ang1/da)+0.5);
               a2 = (int)((ang2/da)+0.5);
               if (r<min_dist){var_energy += forcefield[0][a1][a2]*100*exp(r/min_dist*log(0.01));}
               else{energy += forcefield[dist][a1][a2];}
             }
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
                   dist_n = (r-min_dist)/dr;
                   dist = (int)(dist_n+0.5);
                   a1 = (int)((ang1/da)+0.5);
                   a2 = (int)((ang2/da)+0.5);
                     if (r<min_dist){var_press[numerator] = forcefield[0][a1][a2]*100*exp(r/min_dist*log(0.01));}
                     else{var_press[numerator] = forcefield[dist][a1][a2];}
              numerator++;
              }
              if (press_calc){pressure_X += (var_press[0]-var_press[1])/(delta*4.0)*dist_x;}
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
                    dist_n = (r-min_dist)/dr;
                    dist = (int)(dist_n+0.5);
                    a1 = (int)((ang1/da)+0.5);
                    a2 = (int)((ang2/da)+0.5);
                      if (r<min_dist){var_press[numerator] = forcefield[0][a1][a2]*100*exp(r/min_dist*log(0.01));}
                      else{var_press[numerator] = forcefield[dist][a1][a2];}
               numerator++;
               }
               if (press_calc){pressure_Y += (var_press[0]-var_press[1])/(delta*4.0)*dist_y;}
       }
    }
    en_and_press.energy = energy;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
