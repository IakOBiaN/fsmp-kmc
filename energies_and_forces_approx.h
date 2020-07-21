void charges_coordinates (state &mol, double &Lx, double &Ly)
{
  double d_charges = 4.2; // A
  // for positive charges
  double mol_sin_add_half_carbox_p = mol.sin_phi*0.9659258 + mol.cos_phi*0.258819;
  double mol_cos_add_half_carbox_p = mol.cos_phi*0.9659258 - mol.sin_phi*0.258819;
  double mol_sin_add_half_carbox_p_120 = mol_sin_add_half_carbox_p*(-0.5) + mol_cos_add_half_carbox_p*0.866;
  double mol_cos_add_half_carbox_p_120 = mol_cos_add_half_carbox_p*(-0.5) - mol_sin_add_half_carbox_p*0.866;
  double mol_sin_add_half_carbox_p_240 = mol_sin_add_half_carbox_p*(-0.5) + mol_cos_add_half_carbox_p*(-0.866);
  double mol_cos_add_half_carbox_p_240 = mol_cos_add_half_carbox_p*(-0.5) - mol_sin_add_half_carbox_p*(-0.866);

  // for negative charges
  double mol_sin_add_half_carbox_n = mol.sin_phi*0.9659258 - mol.cos_phi*0.258819;
  double mol_cos_add_half_carbox_n = mol.cos_phi*0.9659258 + mol.sin_phi*0.258819;
  double mol_sin_add_half_carbox_n_120 = mol_sin_add_half_carbox_n*(-0.5) + mol_cos_add_half_carbox_n*0.866;
  double mol_cos_add_half_carbox_n_120 = mol_cos_add_half_carbox_n*(-0.5) - mol_sin_add_half_carbox_n*0.866;
  double mol_sin_add_half_carbox_n_240 = mol_sin_add_half_carbox_n*(-0.5) + mol_cos_add_half_carbox_n*(-0.866);
  double mol_cos_add_half_carbox_n_240 = mol_cos_add_half_carbox_n*(-0.5) - mol_sin_add_half_carbox_n*(-0.866);

  // x,y - coordinates of six charges (3 positive and 3 negative) in the A molecule
  mol.q1x_p = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_p);
  mol.q1y_p = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_p);
  mol.q2x_p = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_p_120);
  mol.q2y_p = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_p_120);
  mol.q3x_p = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_p_240);
  mol.q3y_p = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_p_240);
  mol.q1x_n = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_n);
  mol.q1y_n = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_n);
  mol.q2x_n = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_n_120);
  mol.q2y_n = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_n_120);
  mol.q3x_n = PBC2D(Lx, mol.x + d_charges*mol_cos_add_half_carbox_n_240);
  mol.q3y_n = PBC2D(Ly, mol.y + d_charges*mol_sin_add_half_carbox_n_240);
}

void energy_calculation(state molA, state molB, double &Lx, double &Ly, double &beta, double &r, double &U_LJ, double &U_QQ)
{

	//double derivative_LJ, derivative_QQ;
	double r0 = 7.5877; // Hard core radius in A
	double sigma_r0 = 11.052 - r0; // in A
	double sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
	double eps = 360.0*R; // in J/mol
	const double Ke = 8.987551e+19; // in J*A/C^2
	double N_a = 6.02214e+23;
	double q = 0.45*1.6e-19; // in C
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
	U_LJ += 4*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol

	/// Coulomb's interaction ///
	//derivative_QQ = 0;
	// Interaction energy between 6 charges of the A molecule and 6 charges of the B molecules
	// Positive charges of A molecule and positive charges of B molecule
	dist2 = ((molB.q1x_p - molA.q1x_p)*(molB.q1x_p - molA.q1x_p) + (molB.q1y_p - molA.q1y_p)*(molB.q1y_p - molA.q1y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q1x_p)*(molB.q2x_p - molA.q1x_p) + (molB.q2y_p - molA.q1y_p)*(molB.q2y_p - molA.q1y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q1x_p)*(molB.q3x_p - molA.q1x_p) + (molB.q3y_p - molA.q1y_p)*(molB.q3y_p - molA.q1y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q2x_p)*(molB.q1x_p - molA.q2x_p) + (molB.q1y_p - molA.q2y_p)*(molB.q1y_p - molA.q2y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q2x_p)*(molB.q2x_p - molA.q2x_p) + (molB.q2y_p - molA.q2y_p)*(molB.q2y_p - molA.q2y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q2x_p)*(molB.q3x_p - molA.q2x_p) + (molB.q3y_p - molA.q2y_p)*(molB.q3y_p - molA.q2y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q3x_p)*(molB.q1x_p - molA.q3x_p) + (molB.q1y_p - molA.q3y_p)*(molB.q1y_p - molA.q3y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q3x_p)*(molB.q2x_p - molA.q3x_p) + (molB.q2y_p - molA.q3y_p)*(molB.q2y_p - molA.q3y_p));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q3x_p)*(molB.q3x_p - molA.q3x_p) + (molB.q3y_p - molA.q3y_p)*(molB.q3y_p - molA.q3y_p));
	U_QQ += QQ_const/sqrt(dist2);

	// Negative charges of A molecule and negative charges of B molecule
	dist2 = ((molB.q1x_n - molA.q1x_n)*(molB.q1x_n - molA.q1x_n) + (molB.q1y_n - molA.q1y_n)*(molB.q1y_n - molA.q1y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q1x_n)*(molB.q2x_n - molA.q1x_n) + (molB.q2y_n - molA.q1y_n)*(molB.q2y_n - molA.q1y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q1x_n)*(molB.q3x_n - molA.q1x_n) + (molB.q3y_n - molA.q1y_n)*(molB.q3y_n - molA.q1y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q2x_n)*(molB.q1x_n - molA.q2x_n) + (molB.q1y_n - molA.q2y_n)*(molB.q1y_n - molA.q2y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q2x_n)*(molB.q2x_n - molA.q2x_n) + (molB.q2y_n - molA.q2y_n)*(molB.q2y_n - molA.q2y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q2x_n)*(molB.q3x_n - molA.q2x_n) + (molB.q3y_n - molA.q2y_n)*(molB.q3y_n - molA.q2y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q3x_n)*(molB.q1x_n - molA.q3x_n) + (molB.q1y_n - molA.q3y_n)*(molB.q1y_n - molA.q3y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q3x_n)*(molB.q2x_n - molA.q3x_n) + (molB.q2y_n - molA.q3y_n)*(molB.q2y_n - molA.q3y_n));
	U_QQ += QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q3x_n)*(molB.q3x_n - molA.q3x_n) + (molB.q3y_n - molA.q3y_n)*(molB.q3y_n - molA.q3y_n));
	U_QQ += QQ_const/sqrt(dist2);

	// Negative charges of A molecule and positive charges of B molecule
	dist2 = ((molB.q1x_p - molA.q1x_n)*(molB.q1x_p - molA.q1x_n) + (molB.q1y_p - molA.q1y_n)*(molB.q1y_p - molA.q1y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q1x_n)*(molB.q2x_p - molA.q1x_n) + (molB.q2y_p - molA.q1y_n)*(molB.q2y_p - molA.q1y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q1x_n)*(molB.q3x_p - molA.q1x_n) + (molB.q3y_p - molA.q1y_n)*(molB.q3y_p - molA.q1y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q2x_n)*(molB.q1x_p - molA.q2x_n) + (molB.q1y_p - molA.q2y_n)*(molB.q1y_p - molA.q2y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q2x_n)*(molB.q2x_p - molA.q2x_n) + (molB.q2y_p - molA.q2y_n)*(molB.q2y_p - molA.q2y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q2x_n)*(molB.q3x_p - molA.q2x_n) + (molB.q3y_p - molA.q2y_n)*(molB.q3y_p - molA.q2y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_p - molA.q3x_n)*(molB.q1x_p - molA.q3x_n) + (molB.q1y_p - molA.q3y_n)*(molB.q1y_p - molA.q3y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_p - molA.q3x_n)*(molB.q2x_p - molA.q3x_n) + (molB.q2y_p - molA.q3y_n)*(molB.q2y_p - molA.q3y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_p - molA.q3x_n)*(molB.q3x_p - molA.q3x_n) + (molB.q3y_p - molA.q3y_n)*(molB.q3y_p - molA.q3y_n));
	U_QQ -= QQ_const/sqrt(dist2);

	// Positive charges of A molecule and negative charges of B molecule
	dist2 = ((molB.q1x_n - molA.q1x_p)*(molB.q1x_n - molA.q1x_p) + (molB.q1y_n - molA.q1y_p)*(molB.q1y_n - molA.q1y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q1x_p)*(molB.q2x_n - molA.q1x_p) + (molB.q2y_n - molA.q1y_p)*(molB.q2y_n - molA.q1y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q1x_p)*(molB.q3x_n - molA.q1x_p) + (molB.q3y_n - molA.q1y_p)*(molB.q3y_n - molA.q1y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q2x_p)*(molB.q1x_n - molA.q2x_p) + (molB.q1y_n - molA.q2y_p)*(molB.q1y_n - molA.q2y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q2x_p)*(molB.q2x_n - molA.q2x_p) + (molB.q2y_n - molA.q2y_p)*(molB.q2y_n - molA.q2y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q2x_p)*(molB.q3x_n - molA.q2x_p) + (molB.q3y_n - molA.q2y_p)*(molB.q3y_n - molA.q2y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q1x_n - molA.q3x_p)*(molB.q1x_n - molA.q3x_p) + (molB.q1y_n - molA.q3y_p)*(molB.q1y_n - molA.q3y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q2x_n - molA.q3x_p)*(molB.q2x_n - molA.q3x_p) + (molB.q2y_n - molA.q3y_p)*(molB.q2y_n - molA.q3y_p));
	U_QQ -= QQ_const/sqrt(dist2);

	dist2 = ((molB.q3x_n - molA.q3x_p)*(molB.q3x_n - molA.q3x_p) + (molB.q3y_n - molA.q3y_p)*(molB.q3y_n - molA.q3y_p));
	U_QQ -= QQ_const/sqrt(dist2);
}

results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
		double diff_delta = 0.000001;


		double U_LJ = 0, U_QQ = 0;
		double derivative_LJ, derivative_QQ;
		double r0 = 7.5877; // Hard core radius in A
		double pressure_X = 0, pressure_Y = 0;

		int numerator;

    for (int id = -1; id < 2; id++)
    {
			 dist_x = 0;
			 dist_y = 0;
			 r = 0;
       dist_x = (molB.x - molA.x + id*Lx);
       if (abs(dist_x) > max_dist) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          dist_y = (molB.y - molA.y + jd*Ly);
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
								double t_U_LJ = 0, t_U_QQ = 0;
								energy_calculation(molA, molB, Lx, Ly, beta, r, t_U_LJ, t_U_QQ);

								U_LJ += t_U_LJ;
								U_QQ += t_U_QQ;

								double dist_x_plus_delta, dist_y_plus_delta;
								double t_U_LJ_delta = 0, t_U_QQ_delta = 0;
								state molB_clone = molB;
								molB_clone.x = PBC2D(Lx, molB.x - diff_delta);
								charges_coordinates(molB_clone, Lx, Ly);
								dist_x_plus_delta = molB_clone.x - molA.x + id*Lx;
								r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
								r = sqrt(r2);
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, t_U_LJ_delta, t_U_QQ_delta);
								pressure_X += -((t_U_LJ+t_U_QQ)-(t_U_LJ_delta+t_U_QQ_delta))/diff_delta*dist_x;

								t_U_LJ_delta = 0;
								t_U_QQ_delta = 0;
								molB_clone.x = molB.x;
							 	molB_clone.y = PBC2D(Ly, molB.y - diff_delta);
							  charges_coordinates(molB_clone, Lx, Ly);
								dist_y_plus_delta = molB_clone.y - molA.y + jd*Ly;
								r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
							  r = sqrt(r2);
								energy_calculation(molA, molB_clone, Lx, Ly, beta, r, t_U_LJ_delta, t_U_QQ_delta);
								pressure_Y += -((t_U_LJ+t_U_QQ)-(t_U_LJ_delta+t_U_QQ_delta))/diff_delta*dist_y;
             }
       }
    }
    en_and_press.energy = U_LJ + U_QQ;
		en_and_press.energy_QQ = U_QQ;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
