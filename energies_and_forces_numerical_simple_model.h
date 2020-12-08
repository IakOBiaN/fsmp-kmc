void charges_coordinates (state &mol)
{
  double d_charges = 4.2; // A
  // for positive charges
  double mol_sin_add_half_carbox_p = mol.sin_phi*0.9659258262890682867497431997289 + mol.cos_phi*0.25881904510252076234889883762405;
  double mol_cos_add_half_carbox_p = mol.cos_phi*0.9659258262890682867497431997289 - mol.sin_phi*0.25881904510252076234889883762405;
  double mol_sin_add_half_carbox_p_120 = mol_sin_add_half_carbox_p*(-0.5) + mol_cos_add_half_carbox_p*0.86602540378443864676372317075294;
  double mol_cos_add_half_carbox_p_120 = mol_cos_add_half_carbox_p*(-0.5) - mol_sin_add_half_carbox_p*0.86602540378443864676372317075294;
  double mol_sin_add_half_carbox_p_240 = mol_sin_add_half_carbox_p*(-0.5) + mol_cos_add_half_carbox_p*(-0.86602540378443864676372317075294);
  double mol_cos_add_half_carbox_p_240 = mol_cos_add_half_carbox_p*(-0.5) - mol_sin_add_half_carbox_p*(-0.86602540378443864676372317075294);

  // for negative charges
  double mol_sin_add_half_carbox_n = mol.sin_phi*0.9659258262890682867497431997289 - mol.cos_phi*0.25881904510252076234889883762405;
  double mol_cos_add_half_carbox_n = mol.cos_phi*0.9659258262890682867497431997289 + mol.sin_phi*0.25881904510252076234889883762405;
  double mol_sin_add_half_carbox_n_120 = mol_sin_add_half_carbox_n*(-0.5) + mol_cos_add_half_carbox_n*0.86602540378443864676372317075294;
  double mol_cos_add_half_carbox_n_120 = mol_cos_add_half_carbox_n*(-0.5) - mol_sin_add_half_carbox_n*0.86602540378443864676372317075294;
  double mol_sin_add_half_carbox_n_240 = mol_sin_add_half_carbox_n*(-0.5) + mol_cos_add_half_carbox_n*(-0.86602540378443864676372317075294);
  double mol_cos_add_half_carbox_n_240 = mol_cos_add_half_carbox_n*(-0.5) - mol_sin_add_half_carbox_n*(-0.86602540378443864676372317075294);

  // x,y - coordinates of six charges (3 positive and 3 negative) in the A molecule
	mol.q1x_p = mol.x + d_charges*mol_cos_add_half_carbox_p;
	mol.q2x_p = mol.x + d_charges*mol_cos_add_half_carbox_p_120;
	mol.q3x_p = mol.x + d_charges*mol_cos_add_half_carbox_p_240;
	mol.q1x_n = mol.x + d_charges*mol_cos_add_half_carbox_n;
	mol.q2x_n = mol.x + d_charges*mol_cos_add_half_carbox_n_120;
	mol.q3x_n = mol.x + d_charges*mol_cos_add_half_carbox_n_240;

	mol.q1y_p = mol.y + d_charges*mol_sin_add_half_carbox_p;
  mol.q2y_p = mol.y + d_charges*mol_sin_add_half_carbox_p_120;
  mol.q3y_p = mol.y + d_charges*mol_sin_add_half_carbox_p_240;
  mol.q1y_n = mol.y + d_charges*mol_sin_add_half_carbox_n;
  mol.q2y_n = mol.y + d_charges*mol_sin_add_half_carbox_n_120;
  mol.q3y_n = mol.y + d_charges*mol_sin_add_half_carbox_n_240;
}

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
	else{en = forcefield[dist][a1][a2];}
}

results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
		double diff_delta = 0.000001;


		double r0 = 7.5877; // Hard core radius in A
		state molB_clone = molB;
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
               en_and_press.energy += 1e20;
               continue;
             }
             if (r > min_dist && r <= max_dist)
             {
								double t_U;
								charges_coordinates(molB_clone);
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x, dist_y, t_U);
								energy_cur += t_U;

								// Pressure calculation

								double dist_x_plus_delta, dist_y_plus_delta;
								double t_U_delta = 0;

								molB_clone.x = molB_clone.x - diff_delta;
								dist_x_plus_delta = molB_clone.x - molA.x;
								r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
								r = sqrt(r2);
								charges_coordinates(molB_clone);
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x_plus_delta, dist_y, t_U_delta);
								pressure_X += -(t_U-t_U_delta)/diff_delta*dist_x;
								molB_clone.x = molB_clone.x + diff_delta;

								t_U_delta = 0;
								molB_clone.y = molB_clone.y - diff_delta;
								dist_y_plus_delta = molB_clone.y - molA.y;
								r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
							  r = sqrt(r2);
								charges_coordinates(molB_clone);
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x, dist_y_plus_delta, t_U_delta);
								pressure_Y += -(t_U-t_U_delta)/diff_delta*dist_y;
								molB_clone.y = molB_clone.y + diff_delta;
             }
       }
    }
    en_and_press.energy = energy_cur;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
