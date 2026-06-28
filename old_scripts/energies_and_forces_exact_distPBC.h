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
	mol.q1x_p = mol.x + d_charges*mol_cos_add_half_carbox_p;
  mol.q1y_p = mol.y + d_charges*mol_sin_add_half_carbox_p;
  mol.q2x_p = mol.x + d_charges*mol_cos_add_half_carbox_p_120;
  mol.q2y_p = mol.y + d_charges*mol_sin_add_half_carbox_p_120;
  mol.q3x_p = mol.x + d_charges*mol_cos_add_half_carbox_p_240;
  mol.q3y_p = mol.y + d_charges*mol_sin_add_half_carbox_p_240;
  mol.q1x_n = mol.x + d_charges*mol_cos_add_half_carbox_n;
  mol.q1y_n = mol.y + d_charges*mol_sin_add_half_carbox_n;
  mol.q2x_n = mol.x + d_charges*mol_cos_add_half_carbox_n_120;
  mol.q2y_n = mol.y + d_charges*mol_sin_add_half_carbox_n_120;
  mol.q3x_n = mol.x + d_charges*mol_cos_add_half_carbox_n_240;
  mol.q3y_n = mol.y + d_charges*mol_sin_add_half_carbox_n_240;;
}

results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
  	double ang1,ang2;
  	double r;	// Distance between A and B molecules
  	double r2;	// Distance sqaured
  	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies


		double U_LJ = 0, U_QQ = 0;
		double derivative_LJ, derivative_QQ_x, derivative_QQ_y;
		double r0 = 7.5877; // Hard core radius in A
		double sigma_r0 = 11.052 - r0; // in A
		double eps = 360.0*R; // in J/mol
		const double Ke = 8.987551e+19; // in J*A/C^2
		double N_a = 6.02214e+23;
		double q = 0.45*1.6e-19; // in C
		double q2 = q*q;
		double QQ_const = Ke*q2*N_a;
		double pressure_X = 0, pressure_Y = 0;
		double dist2 = 0;

		/// Lennard-Johns interaction ///
		double r_r0 = 0;
		double sigma_r0_6 = 0;
		double r_r0_6 = 0;
		double lj_dist_6 = 0;

		dist_x = 0;
		dist_y = 0;
		r = 0;

		dist_x = distPBC(Lx, (molB.x - molA.x));
		if (abs(dist_x) < max_dist)
		{
			dist_y = distPBC(Ly, (molB.y - molA.y));
			if (abs(dist_y) < max_dist)
			{
				r2 = dist_x*dist_x + dist_y*dist_y;
				if (r2 <= min_dist*min_dist)
					{
						en_and_press.energy += 1e20;
						return en_and_press;
					}
				if (r2 > min_dist*min_dist && r2 <= max_dist*max_dist)
					{
						r = sqrt(r2);
						derivative_LJ = 0;
						/// Lennard-Johns interaction ///
						r_r0 = r - r0;
						sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
						r_r0_6 = r_r0*r_r0*r_r0*r_r0*r_r0*r_r0;
						lj_dist_6 = sigma_r0_6/r_r0_6;
						U_LJ += 4*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol
						derivative_LJ += 24*eps*lj_dist_6 * (1-2.0*lj_dist_6)/r_r0;

						/// Coulomb's interaction ///
						// Interaction energy between 6 charges of the A molecule and 6 charges of the B molecules
						double dist_ab_x = 0;
						double dist_ab_y = 0;
						double u_columb = 0;
						derivative_QQ_x = 0;
						derivative_QQ_y = 0;

						// Positive charges of A molecule and positive charges of B molecule
						dist_ab_x = (molB.q1x_p - molA.q1x_p);
						dist_ab_y = (molB.q1y_p - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q1x_p);
						dist_ab_y = (molB.q2y_p - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q1x_p);
						dist_ab_y = (molB.q3y_p - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_p - molA.q2x_p);
						dist_ab_y = (molB.q1y_p - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q2x_p);
						dist_ab_y = (molB.q2y_p - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q2x_p);
						dist_ab_y = (molB.q3y_p - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_p - molA.q3x_p);
						dist_ab_y = (molB.q1y_p - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q3x_p);
						dist_ab_y = (molB.q2y_p - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q3x_p);
						dist_ab_y = (molB.q3y_p - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						// Negative charges of A molecule and negative charges of B molecule
						dist_ab_x = (molB.q1x_n - molA.q1x_n);
						dist_ab_y = (molB.q1y_n - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_n - molA.q1x_n);
						dist_ab_y = (molB.q2y_n - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q1x_n);
						dist_ab_y = (molB.q3y_n - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_n - molA.q2x_n);
						dist_ab_y = (molB.q1y_n - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;;

						dist_ab_x = (molB.q2x_n - molA.q2x_n);
						dist_ab_y = (molB.q2y_n - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q2x_n);
						dist_ab_y = (molB.q3y_n - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_n - molA.q3x_n);
						dist_ab_y = (molB.q1y_n - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_n - molA.q3x_n);
						dist_ab_y = (molB.q2y_n - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q3x_n);
						dist_ab_y = (molB.q3y_n - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						// Negative charges of A molecule and positive charges of B molecule
						dist_ab_x = (molB.q1x_p - molA.q1x_n);
						dist_ab_y = (molB.q1y_p - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q1x_n);
						dist_ab_y = (molB.q2y_p - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q1x_n);
						dist_ab_y = (molB.q3y_p - molA.q1y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_p - molA.q2x_n);
						dist_ab_y = (molB.q1y_p - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q2x_n);
						dist_ab_y = (molB.q2y_p - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q2x_n);
						dist_ab_y = (molB.q3y_p - molA.q2y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_p - molA.q3x_n);
						dist_ab_y = (molB.q1y_p - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_p - molA.q3x_n);
						dist_ab_y = (molB.q2y_p - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_p - molA.q3x_n);
						dist_ab_y = (molB.q3y_p - molA.q3y_n);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						// Positive charges of A molecule and negative charges of B molecule
						dist_ab_x = (molB.q1x_n - molA.q1x_p);
						dist_ab_y = (molB.q1y_n - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_n - molA.q1x_p);
						dist_ab_y = (molB.q2y_n - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q1x_p);
						dist_ab_y = (molB.q3y_n - molA.q1y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_n - molA.q2x_p);
						dist_ab_y = (molB.q1y_n - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_n - molA.q2x_p);
						dist_ab_y = (molB.q2y_n - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q2x_p);
						dist_ab_y = (molB.q3y_n - molA.q2y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q1x_n - molA.q3x_p);
						dist_ab_y = (molB.q1y_n - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q2x_n - molA.q3x_p);
						dist_ab_y = (molB.q2y_n - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						dist_ab_x = (molB.q3x_n - molA.q3x_p);
						dist_ab_y = (molB.q3y_n - molA.q3y_p);
						dist_ab_x = distPBC(Lx, dist_ab_x);
						dist_ab_y = distPBC(Ly, dist_ab_y);
						dist2 = (dist_ab_x*dist_ab_x + dist_ab_y*dist_ab_y);
						u_columb = -QQ_const/sqrt(dist2);
						derivative_QQ_x += -u_columb*dist_ab_x/dist2;
						derivative_QQ_y += -u_columb*dist_ab_y/dist2;
						U_QQ += u_columb;

						pressure_X += -derivative_LJ*dist_x*dist_x/r - derivative_QQ_x*dist_x;
						pressure_Y += -derivative_LJ*dist_y*dist_y/r - derivative_QQ_y*dist_y;
					}
				}
       }

    en_and_press.energy = U_LJ + U_QQ;
		en_and_press.energy_QQ = U_QQ;
    en_and_press.p_X = pressure_X;
    en_and_press.p_Y = pressure_Y;
    return en_and_press;
}
