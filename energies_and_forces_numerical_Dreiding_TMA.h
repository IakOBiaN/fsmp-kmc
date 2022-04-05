void energy_calculation(state molA, state molB, double &Lx, double &Ly, double &beta, double &r, double dist_x, double dist_y, double &en)
{

	double ang_molA = molA.phi;
	double ang_molB = molB.phi;
	double dist_n; // Float index in the numerical potential
	int dist,a1,a2; // indexes in the numerical potential array
	double ang1,ang2;
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
	if (r<min_dist){en = 1e20;}
	else{en = forcefield[dist][a1][a2]*molA.damping_coeff*molB.damping_coeff+molA.ex_field_coeff+molB.ex_field_coeff;}
}

results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta, bool pressure_calc)
{
    results en_and_press;
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
		double diff_delta = 0.01;


//		double r0 = 7.5877; // Hard core radius in A
		state molB_clone = molB;
		state molA_clone = molA;
		double energy_cur = 0;
		double pressure_X = 0, pressure_Y = 0;

    for (int id = -1; id < 2; id++)
    {
			 dist_x = 0;
			 dist_y = 0;
			 r = 0;
			 molB_clone.x = molB.x + id*Lx;
       dist_x = molB_clone.x - molA.x;
       if (abs(dist_x) > max_dist) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
				  molB_clone.y = molB.y + jd*Ly;
          dist_y = molB_clone.y - molA.y;
          if (abs(dist_y) > max_dist) {continue;}
             r2 = dist_x*dist_x + dist_y*dist_y;
             r = sqrt(r2);
             if (r <= min_dist)
             {
               HC_radius = true;
               en_and_press.energy += 1e20;
               continue;
             }
             if (r > min_dist && r <= (max_dist-2.0*diff_delta))
             {
								double t_U;
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x, dist_y, t_U);
								energy_cur += t_U;
								if (pressure_calc)
								{
									// Pressure calculation
									double dist_x_plus_delta, dist_y_plus_delta;
									double t_U_delta_1 = 0, t_U_delta_2 = 0;

									molB_clone.x = molB_clone.x - diff_delta;
									molA_clone.x = molA_clone.x + diff_delta;
									dist_x_plus_delta = molB_clone.x - molA_clone.x;
									r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
									r = sqrt(r2);
									energy_calculation(molA_clone, molB_clone, Lx, Ly, beta, r, dist_x_plus_delta, dist_y, t_U_delta_1);

									molB_clone.x = molB_clone.x + 2.0*diff_delta;
									molA_clone.x = molA_clone.x - 2.0*diff_delta;
									dist_x_plus_delta = molB_clone.x - molA_clone.x;
									r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
									r = sqrt(r2);
									energy_calculation(molA_clone, molB_clone, Lx, Ly, beta, r, dist_x_plus_delta, dist_y, t_U_delta_2);

									pressure_X += -(t_U_delta_2-t_U_delta_1)/(4.0*diff_delta)*dist_x;
									molB_clone.x = molB_clone.x - diff_delta;
									molA_clone.x = molA_clone.x + diff_delta;

									t_U_delta_1 = 0;
									t_U_delta_2 = 0;

									molB_clone.y = molB_clone.y - diff_delta;
									molA_clone.y = molA_clone.y + diff_delta;
									dist_y_plus_delta = molB_clone.y - molA_clone.y;
									r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
								  r = sqrt(r2);
									energy_calculation(molA_clone, molB_clone, Lx, Ly, beta, r, dist_x, dist_y_plus_delta, t_U_delta_1);

									molB_clone.y = molB_clone.y + 2.0*diff_delta;
									molA_clone.y = molA_clone.y - 2.0*diff_delta;
									dist_y_plus_delta = molB_clone.y - molA_clone.y;
									r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
								  r = sqrt(r2);
									energy_calculation(molA_clone, molB_clone, Lx, Ly, beta, r, dist_x, dist_y_plus_delta, t_U_delta_2);

									pressure_Y += -(t_U_delta_2-t_U_delta_1)/(4.0*diff_delta)*dist_y;
									molB_clone.y = molB_clone.y - diff_delta;
									molA_clone.y = molA_clone.y + diff_delta;
								}
             }
       }
    }
    en_and_press.energy = energy_cur;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
