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

void energy_calculation(state molA, state molB, double &Lx, double &Ly, double &beta, double &r, double &t_U_LJ, double &t_U_QQ)
{

	//double derivative_LJ, derivative_QQ;
	double r0 = 7.5877; // Hard core radius in A
	double sigma_r0 = 11.052 - r0; // in A
	double sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
	double eps = 360.0*R; // in J/mol
	const double Ke = 8.9875517923e+19; // in J*A/C^2
	double N_a = 6.02214076e+23;
	double q = 0.45*1.602176634e-19; // in C
	double q2 = q*q;
	double QQ_const = Ke*q2*N_a;
	double dist2 = 0;

	/// Lennard-Johns interaction ///
	double r_r0 = 0;
	double r_r0_6 = 0;
	double lj_dist_6 = 0;

	//derivative_LJ = 0;
	/// Lennard-Johns interaction ///
	r_r0 = r - r0;
	r_r0_6 = r_r0*r_r0*r_r0*r_r0*r_r0*r_r0;
	lj_dist_6 = sigma_r0_6/r_r0_6;
	t_U_LJ += 4*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol

	// Damping field effect
	t_U_LJ *= molA.damping_coeff*molB.damping_coeff;

	/// Coulomb's interaction ///
	//derivative_QQ = 0;
	// Interaction energy between 6 charges of the A molecule and 6 charges of the B molecules
	// Positive charges of A molecule and positive charges of B molecule
	dist2 = ((molB.q1x_p - molA.q1x_p)*(molB.q1x_p - molA.q1x_p) + (molB.q1y_p - molA.q1y_p)*(molB.q1y_p - molA.q1y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q1x_p)*(molB.q2x_p - molA.q1x_p) + (molB.q2y_p - molA.q1y_p)*(molB.q2y_p - molA.q1y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q1x_p)*(molB.q3x_p - molA.q1x_p) + (molB.q3y_p - molA.q1y_p)*(molB.q3y_p - molA.q1y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q2x_p)*(molB.q1x_p - molA.q2x_p) + (molB.q1y_p - molA.q2y_p)*(molB.q1y_p - molA.q2y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q2x_p)*(molB.q2x_p - molA.q2x_p) + (molB.q2y_p - molA.q2y_p)*(molB.q2y_p - molA.q2y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q2x_p)*(molB.q3x_p - molA.q2x_p) + (molB.q3y_p - molA.q2y_p)*(molB.q3y_p - molA.q2y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q3x_p)*(molB.q1x_p - molA.q3x_p) + (molB.q1y_p - molA.q3y_p)*(molB.q1y_p - molA.q3y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q3x_p)*(molB.q2x_p - molA.q3x_p) + (molB.q2y_p - molA.q3y_p)*(molB.q2y_p - molA.q3y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q3x_p)*(molB.q3x_p - molA.q3x_p) + (molB.q3y_p - molA.q3y_p)*(molB.q3y_p - molA.q3y_p));
	t_U_QQ += QQ_const/sqrt(dist2);

	// Negative charges of A molecule and negative charges of B molecule
	dist2 = ((molB.q1x_n - molA.q1x_n)*(molB.q1x_n - molA.q1x_n) + (molB.q1y_n - molA.q1y_n)*(molB.q1y_n - molA.q1y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q1x_n)*(molB.q2x_n - molA.q1x_n) + (molB.q2y_n - molA.q1y_n)*(molB.q2y_n - molA.q1y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q1x_n)*(molB.q3x_n - molA.q1x_n) + (molB.q3y_n - molA.q1y_n)*(molB.q3y_n - molA.q1y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q2x_n)*(molB.q1x_n - molA.q2x_n) + (molB.q1y_n - molA.q2y_n)*(molB.q1y_n - molA.q2y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q2x_n)*(molB.q2x_n - molA.q2x_n) + (molB.q2y_n - molA.q2y_n)*(molB.q2y_n - molA.q2y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q2x_n)*(molB.q3x_n - molA.q2x_n) + (molB.q3y_n - molA.q2y_n)*(molB.q3y_n - molA.q2y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q3x_n)*(molB.q1x_n - molA.q3x_n) + (molB.q1y_n - molA.q3y_n)*(molB.q1y_n - molA.q3y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q3x_n)*(molB.q2x_n - molA.q3x_n) + (molB.q2y_n - molA.q3y_n)*(molB.q2y_n - molA.q3y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q3x_n)*(molB.q3x_n - molA.q3x_n) + (molB.q3y_n - molA.q3y_n)*(molB.q3y_n - molA.q3y_n));
	t_U_QQ += QQ_const/sqrt(dist2);

	// Negative charges of A molecule and positive charges of B molecule
	dist2 = ((molB.q1x_p - molA.q1x_n)*(molB.q1x_p - molA.q1x_n) + (molB.q1y_p - molA.q1y_n)*(molB.q1y_p - molA.q1y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q1x_n)*(molB.q2x_p - molA.q1x_n) + (molB.q2y_p - molA.q1y_n)*(molB.q2y_p - molA.q1y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q1x_n)*(molB.q3x_p - molA.q1x_n) + (molB.q3y_p - molA.q1y_n)*(molB.q3y_p - molA.q1y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q2x_n)*(molB.q1x_p - molA.q2x_n) + (molB.q1y_p - molA.q2y_n)*(molB.q1y_p - molA.q2y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q2x_n)*(molB.q2x_p - molA.q2x_n) + (molB.q2y_p - molA.q2y_n)*(molB.q2y_p - molA.q2y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q2x_n)*(molB.q3x_p - molA.q2x_n) + (molB.q3y_p - molA.q2y_n)*(molB.q3y_p - molA.q2y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q3x_n)*(molB.q1x_p - molA.q3x_n) + (molB.q1y_p - molA.q3y_n)*(molB.q1y_p - molA.q3y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q3x_n)*(molB.q2x_p - molA.q3x_n) + (molB.q2y_p - molA.q3y_n)*(molB.q2y_p - molA.q3y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q3x_n)*(molB.q3x_p - molA.q3x_n) + (molB.q3y_p - molA.q3y_n)*(molB.q3y_p - molA.q3y_n));
	t_U_QQ -= QQ_const/sqrt(dist2);

	// Positive charges of A molecule and negative charges of B molecule
	dist2 = ((molB.q1x_n - molA.q1x_p)*(molB.q1x_n - molA.q1x_p) + (molB.q1y_n - molA.q1y_p)*(molB.q1y_n - molA.q1y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q1x_p)*(molB.q2x_n - molA.q1x_p) + (molB.q2y_n - molA.q1y_p)*(molB.q2y_n - molA.q1y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q1x_p)*(molB.q3x_n - molA.q1x_p) + (molB.q3y_n - molA.q1y_p)*(molB.q3y_n - molA.q1y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q2x_p)*(molB.q1x_n - molA.q2x_p) + (molB.q1y_n - molA.q2y_p)*(molB.q1y_n - molA.q2y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q2x_p)*(molB.q2x_n - molA.q2x_p) + (molB.q2y_n - molA.q2y_p)*(molB.q2y_n - molA.q2y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q2x_p)*(molB.q3x_n - molA.q2x_p) + (molB.q3y_n - molA.q2y_p)*(molB.q3y_n - molA.q2y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q3x_p)*(molB.q1x_n - molA.q3x_p) + (molB.q1y_n - molA.q3y_p)*(molB.q1y_n - molA.q3y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q3x_p)*(molB.q2x_n - molA.q3x_p) + (molB.q2y_n - molA.q3y_p)*(molB.q2y_n - molA.q3y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q3x_p)*(molB.q3x_n - molA.q3x_p) + (molB.q3y_n - molA.q3y_p)*(molB.q3y_n - molA.q3y_p));
	t_U_QQ -= QQ_const/sqrt(dist2);

	t_U_QQ *= molA.damping_coeff*molB.damping_coeff;
}

results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta, bool abcd)
{
    results en_and_press;
    if ((molA.damping_coeff == 0) || (molB.damping_coeff == 0))
    {
      return en_and_press;
    }
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
		double dist_x_2, dist_y_2;
		double diff_delta = 0.02;


		double U_LJ = 0, U_QQ = 0;
		double pressure_X = 0, pressure_Y = 0;
		state molB_clone = molB;

		double dist_x_plus_delta, dist_y_plus_delta;
		double dist_x_plus_delta_2, dist_y_plus_delta_2;
		double t_U_LJ_delta, t_U_QQ_delta;

		for (int id = -1; id < 2; id++)
//		for (int id = 0; id < 1; id++)
    {
			 if (HC_radius) {break;}
			 molB_clone.x = molB.x + id*Lx;
       dist_x = molB_clone.x - molA.x;
			 dist_x_2 = dist_x*dist_x;
       if (dist_x_2 > max_dist_2) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
				  if (HC_radius) {break;}
				  molB_clone.y = molB.y + jd*Ly;
          dist_y = molB_clone.y - molA.y;
					dist_y_2 = dist_y*dist_y;
          if (dist_y_2 > max_dist_2) {continue;}
          r2 = dist_x_2 + dist_y_2;
             if (r2 <= min_dist_2)
             {
               en_and_press.energy += E_INF/beta;
               HC_radius = true;
               break;
             }
					if (r2 > min_dist_2 && r2 <= max_dist_2)
					{
					 	r = sqrt(r2);
						double t_U_LJ = 0, t_U_QQ = 0;
						charges_coordinates(molB_clone);
						energy_calculation(molA, molB_clone, Lx, Ly, beta, r, t_U_LJ, t_U_QQ);
						U_LJ += t_U_LJ;
						U_QQ += t_U_QQ;

						t_U_LJ_delta = 0, t_U_QQ_delta = 0;
						molB_clone.x = molB_clone.x - diff_delta;
						dist_x_plus_delta = molB_clone.x - molA.x;
						r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
						r = sqrt(r2);
						charges_coordinates(molB_clone);
						molB_clone.damping_coeff = damping_field(molB_clone.x, Lx);
						molB_clone.ex_field_coeff = external_field(molB_clone.x, Lx);
						energy_calculation(molA, molB_clone, Lx, Ly, beta, r, t_U_LJ_delta, t_U_QQ_delta);
						pressure_X += -((t_U_LJ+t_U_QQ)-(t_U_LJ_delta+t_U_QQ_delta))/diff_delta*dist_x;
						molB_clone.x = molB_clone.x + diff_delta;
						charges_coordinates(molB_clone);
						molB_clone.damping_coeff = damping_field(molB_clone.x, Lx);
						molB_clone.ex_field_coeff = external_field(molB_clone.x, Lx);

						t_U_LJ_delta = 0;
						t_U_QQ_delta = 0;
						molB_clone.y = molB_clone.y - diff_delta;
						dist_y_plus_delta = molB_clone.y - molA.y;
						r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
						r = sqrt(r2);
						charges_coordinates(molB_clone);
						energy_calculation(molA, molB_clone, Lx, Ly, beta, r, t_U_LJ_delta, t_U_QQ_delta);
						pressure_Y += -((t_U_LJ+t_U_QQ)-(t_U_LJ_delta+t_U_QQ_delta))/diff_delta*dist_y;
						molB_clone.y = molB_clone.y + diff_delta;
					}
       }
    }
		if (HC_radius == false)
		{
			en_and_press.energy = U_LJ + U_QQ;
			en_and_press.p_X = pressure_X;
			en_and_press.p_Y = pressure_Y;
	    if(en_and_press.energy >= E_INF/beta)
				{
					en_and_press.energy = E_INF/beta;
					en_and_press.p_X = 0;
					en_and_press.p_Y = 0;
				}
		}
		else {en_and_press.energy = E_INF/beta; en_and_press.p_X = 0; en_and_press.p_Y = 0;}
    return en_and_press;

}


int check_HC(state molA, state molB, double &Lx, double &Ly)
{
      if ((molA.damping_coeff == 0) || (molB.damping_coeff == 0))
      {
        return 0;
      }
	  	double r;	// Distance between A and B molecules
	  	double r2;	// Distance sqaured
	  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
			double dist_x_2, dist_y_2;

	//		double r0 = 7.5877; // Hard core radius in A
			state molB_clone = molB;

	    for (int id = -1; id < 2; id++)
	    {
				 dist_x = 0;
				 dist_y = 0;
				 r = 0;
				 molB_clone.x = molB.x + id*Lx;
	       dist_x = molB_clone.x - molA.x;
				 dist_x_2 = dist_x*dist_x;
	       if (dist_x_2 > max_dist_2) {continue;}
	       for (int jd = -1; jd < 2; jd++)
	       {
					  molB_clone.y = molB.y + jd*Ly;
	          dist_y = molB_clone.y - molA.y;
						dist_y_2 = dist_y*dist_y;
	          if (dist_y_2 > max_dist_2) {continue;}
	          r2 = dist_x_2 + dist_y_2;
	             if (r2 <= min_dist_2)
	             {
	               HC_radius = true;
	               break;
	             }
	       }
	    }

      return 0;
}
