results coarse_grain_potential(state molA, state molB, double &Lx, double &Ly, bool pressure_calculation)
{
	double max_dist = 50.0;
	double alpha = 15.0; // Half of angle between oxygens in the carboxylic group of TMA molecule
	double gamma = 120.0;
	double double_gamma = 240.0;
	double r0 = 7.5877; // Hard core radius in A
	double sigma_r0 = 11.052 - 7.5877; // in A
	double eps = 360.0*R; // in J/mol
	const double Ke = 8.987551e+19; // in J*A/C^2
	double N_a = 6.02214e+23;
	double q = 0.45*1.6e-19; // in C
	double q2 = q*q;
	double d_charges = 4.2; // A
	double QQ_const = Ke*q2*N_a;
	results en_and_press;
	double ang1,ang2;
	double r;	// Distance between A and B molecules
	double r2;	// Distance sqaured
	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
	double energy = 0, var_energy = 0, pressure_X = 0, pressure_Y = 0;
	double dist_x_plus_delta, dist_y_plus_delta; // Distance between A and B molecules including derivation step
	double molA_x  = molA.x;
	double molA_y  = molA.y;
	double molB_x = molB.x;
	double molB_y = molB.y;
	double ang_molA = molA.phi;
	double ang_molB = molB.phi;
	double sin_ang_molA = molA.sin_phi;
	double cos_ang_molA = molA.cos_phi;
	double sin_ang_molB = molB.sin_phi;
	double cos_ang_molB = molB.cos_phi;
	double U_LJ=0;
	double U_QQ=0;

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
						double dang = dist_x/r;	// Calculate the cosine of the angle between OX and distance vector
						if (dist_y<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
						ang1 = ang_molA - dang; // Angle of A molecule when distance vector is parallel to the OX axis
						ang2 = ang_molB - dang; // Angle of A molecule when distance vector is parallel to the OX axis
						// Molecule should always has the angle in the range of 0-360 degrees
						if (ang1<0) {ang1 += 360.0;}
						if (ang2<0) {ang2 += 360.0;}
						if (ang1>359.999) {ang1 -= 360.0;}
						if (ang2>359.999) {ang2 -= 360.0;}

						/// Lennard-Johns interaction ///
						double r_r0 = r - r0;
						double sigma_r0_6 = sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0*sigma_r0;
						double r_r0_6 = r_r0*r_r0*r_r0*r_r0*r_r0*r_r0;
						double lj_dist_6 = sigma_r0_6/r_r0_6;
						U_LJ += 4.0*eps*(lj_dist_6 * (lj_dist_6 - 1.0)); // in J/mol

						/// Coulomb's interaction ///
						double sin_half_carbox = sin(alpha/180.0*PI);
						double cos_half_carbox = cos(alpha/180.0*PI);

						// Coordinates of charges in the molecule A
						double sin_ang1 = sin(ang1/180.0*PI);
						double cos_ang1 = cos(ang1/180.0*PI);
						// for positive charges
						double molA_sin_add_half_carbox_p = sin_ang1*cos_half_carbox + cos_ang1*sin_half_carbox;
						double molA_cos_add_half_carbox_p = cos_ang1*cos_half_carbox - sin_ang1*sin_half_carbox;
						double molA_sin_add_half_carbox_p_gamma = molA_sin_add_half_carbox_p*cos(gamma/180.0*PI) + molA_cos_add_half_carbox_p*sin(gamma/180.0*PI);
						double molA_cos_add_half_carbox_p_gamma = molA_cos_add_half_carbox_p*cos(gamma/180.0*PI) - molA_sin_add_half_carbox_p*sin(gamma/180.0*PI);
						double molA_sin_add_half_carbox_p_double_gamma = molA_sin_add_half_carbox_p*cos(double_gamma/180.0*PI) + molA_cos_add_half_carbox_p*sin(double_gamma/180.0*PI);
						double molA_cos_add_half_carbox_p_double_gamma = molA_cos_add_half_carbox_p*cos(double_gamma/180.0*PI) - molA_sin_add_half_carbox_p*sin(double_gamma/180.0*PI);

						// for negative charges
						double molA_sin_add_half_carbox_n = sin_ang1*cos_half_carbox - cos_ang1*sin_half_carbox;
						double molA_cos_add_half_carbox_n = cos_ang1*cos_half_carbox + sin_ang1*sin_half_carbox;
						double molA_sin_add_half_carbox_n_gamma = molA_sin_add_half_carbox_n*cos(gamma/180.0*PI) + molA_cos_add_half_carbox_n*sin(gamma/180.0*PI);
						double molA_cos_add_half_carbox_n_gamma = molA_cos_add_half_carbox_n*cos(gamma/180.0*PI) - molA_sin_add_half_carbox_n*sin(gamma/180.0*PI);
						double molA_sin_add_half_carbox_n_double_gamma = molA_sin_add_half_carbox_n*cos(double_gamma/180.0*PI) + molA_cos_add_half_carbox_n*sin(double_gamma/180.0*PI);
						double molA_cos_add_half_carbox_n_double_gamma = molA_cos_add_half_carbox_n*cos(double_gamma/180.0*PI) - molA_sin_add_half_carbox_n*sin(double_gamma/180.0*PI);

						// x,y - coordinates of six charges (3 positive and 3 negative) in the A molecule
						double molA_x_q1_p = molA_x + d_charges*molA_cos_add_half_carbox_p;
						double molA_y_q1_p = molA_y + d_charges*molA_sin_add_half_carbox_p;
						double molA_x_q2_p = molA_x + d_charges*molA_cos_add_half_carbox_p_gamma;
						double molA_y_q2_p = molA_y + d_charges*molA_sin_add_half_carbox_p_gamma;
						double molA_x_q3_p = molA_x + d_charges*molA_cos_add_half_carbox_p_double_gamma;
						double molA_y_q3_p = molA_y + d_charges*molA_sin_add_half_carbox_p_double_gamma;
						double molA_x_q1_n = molA_x + d_charges*molA_cos_add_half_carbox_n;
						double molA_y_q1_n = molA_y + d_charges*molA_sin_add_half_carbox_n;
						double molA_x_q2_n = molA_x + d_charges*molA_cos_add_half_carbox_n_gamma;
						double molA_y_q2_n = molA_y + d_charges*molA_sin_add_half_carbox_n_gamma;
						double molA_x_q3_n = molA_x + d_charges*molA_cos_add_half_carbox_n_double_gamma;
						double molA_y_q3_n = molA_y + d_charges*molA_sin_add_half_carbox_n_double_gamma;

						// Coordinates of charges in the molecule B
						double sin_ang2 = sin(ang2/180.0*PI);
						double cos_ang2 = cos(ang2/180.0*PI);
						// for positive charges
						double molB_sin_add_half_carbox_p = sin_ang2*cos_half_carbox + cos_ang2*sin_half_carbox;
						double molB_cos_add_half_carbox_p = cos_ang2*cos_half_carbox - sin_ang2*sin_half_carbox;
						double molB_sin_add_half_carbox_p_gamma = molB_sin_add_half_carbox_p*cos(gamma/180.0*PI) + molB_cos_add_half_carbox_p*sin(gamma/180.0*PI);
						double molB_cos_add_half_carbox_p_gamma = molB_cos_add_half_carbox_p*cos(gamma/180.0*PI) - molB_sin_add_half_carbox_p*sin(gamma/180.0*PI);
						double molB_sin_add_half_carbox_p_double_gamma = molB_sin_add_half_carbox_p*cos(double_gamma/180.0*PI) + molB_cos_add_half_carbox_p*sin(double_gamma/180.0*PI);
						double molB_cos_add_half_carbox_p_double_gamma = molB_cos_add_half_carbox_p*cos(double_gamma/180.0*PI) - molB_sin_add_half_carbox_p*sin(double_gamma/180.0*PI);
						// for negative charges
						double molB_sin_add_half_carbox_n = sin_ang2*cos_half_carbox - cos_ang2*sin_half_carbox;
						double molB_cos_add_half_carbox_n = cos_ang2*cos_half_carbox + sin_ang2*sin_half_carbox;
						double molB_sin_add_half_carbox_n_gamma = molB_sin_add_half_carbox_n*cos(gamma/180.0*PI) + molB_cos_add_half_carbox_n*sin(gamma/180.0*PI);
						double molB_cos_add_half_carbox_n_gamma = molB_cos_add_half_carbox_n*cos(gamma/180.0*PI) - molB_sin_add_half_carbox_n*sin(gamma/180.0*PI);
						double molB_sin_add_half_carbox_n_double_gamma = molB_sin_add_half_carbox_n*cos(double_gamma/180.0*PI) + molB_cos_add_half_carbox_n*sin(double_gamma/180.0*PI);
						double molB_cos_add_half_carbox_n_double_gamma = molB_cos_add_half_carbox_n*cos(double_gamma/180.0*PI) - molB_sin_add_half_carbox_n*sin(double_gamma/180.0*PI);
						// x,y - coordinates of six charges (3 positive and 3 negative) in the B molecule
						double molB_x_q1_p = molB_x + d_charges*molB_cos_add_half_carbox_p;
						double molB_y_q1_p = molB_y + d_charges*molB_sin_add_half_carbox_p;
						double molB_x_q2_p = molB_x + d_charges*molB_cos_add_half_carbox_p_gamma;
						double molB_y_q2_p = molB_y + d_charges*molB_sin_add_half_carbox_p_gamma;
						double molB_x_q3_p = molB_x + d_charges*molB_cos_add_half_carbox_p_double_gamma;
						double molB_y_q3_p = molB_y + d_charges*molB_sin_add_half_carbox_p_double_gamma;
						double molB_x_q1_n = molB_x + d_charges*molB_cos_add_half_carbox_n;
						double molB_y_q1_n = molB_y + d_charges*molB_sin_add_half_carbox_n;
						double molB_x_q2_n = molB_x + d_charges*molB_cos_add_half_carbox_n_gamma;
						double molB_y_q2_n = molB_y + d_charges*molB_sin_add_half_carbox_n_gamma;
						double molB_x_q3_n = molB_x + d_charges*molB_cos_add_half_carbox_n_double_gamma;
						double molB_y_q3_n = molB_y + d_charges*molB_sin_add_half_carbox_n_double_gamma;

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

						if (pressure_calculation)
						{
						}
					}
		}
	}
//	U_QQ = 0;
	U_LJ = 0;
	en_and_press.energy = U_QQ + U_LJ;
//	en_and_press.p_X = pressure_X;
//	en_and_press.p_Y = pressure_Y;
	return en_and_press;

}
