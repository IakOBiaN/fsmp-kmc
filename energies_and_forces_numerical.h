void energy_calculation(const state &molA, const state &molB, double & /*Lx*/, double & /*Ly*/, double &beta, double &r, double dist_x, double dist_y, double &en)
{
	bool tail_correction = true; // Use or not the tail correction
	bool interpol = true; // Use or not the interpolation procedure
	double ang_molA = molA.phi;
	double ang_molB = molB.phi;
	double dist_n; // Float index in the numerical potential
	int dist,a1,a2; // indexes in the numerical potential array
	double ang1,ang2;
	double dang = dist_x/r;	// Calculate the cosine of the angle between OX and distance vector

	double u_cut = 0;
	double u_before_cut = 0;

	if (dist_y<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
	ang1 = ang_molA - dang;
	ang2 = ang_molB - dang;
	// Reduce the orientation into the stored angular range. For an unfolded grid that is
	// [0, 360); for a grid folded to the molecule's symmetry period it is [0, ff_fold_deg),
	// which is equivalent by that rotational symmetry.
	if (ff_fold_deg >= 359.999)
	{
		if (ang1 < 0) {ang1 += 360.0;}
		if (ang2 < 0) {ang2 += 360.0;}
		if (ang1 > 359.95) {ang1 -= 360.0;}
		if (ang2 > 359.95) {ang2 -= 360.0;}
	}
	else
	{
		ang1 = fmod(ang1, ff_fold_deg); if (ang1 < 0) {ang1 += ff_fold_deg;}
		ang2 = fmod(ang2, ff_fold_deg); if (ang2 < 0) {ang2 += ff_fold_deg;}
	}
	dist_n = (r - min_dist) / dr;
	dist = (int)(dist_n + 0.5);
	a1 = (int)((ang1 / da) + 0.5);
	a2 = (int)((ang2 / da) + 0.5);
	if(r <= min_dist) {en = (E_INF / beta) * molA.damping_coeff * molB.damping_coeff;}
	else
			{
				if (r > max_dist) {en = 0;}
					else
							{
								if (tail_correction)
									{
										u_cut = FF(cut_index, a1, a2);
										u_before_cut = FF(cut_index - 1, a1, a2);
										if (interpol)
										{
											en = interpolation(dist_n, ang1, ang2, beta) - u_cut + (u_cut - u_before_cut) * (r - max_dist) / dr;
										}
										else
										{
											en = FF(dist, a1, a2) - u_cut + (u_cut - u_before_cut) * (r - max_dist) / dr;
										}
									}
								else
									{
										if (interpol)
										{
											en = interpolation(dist_n, ang1, ang2, beta);
										}
										else
										{
											en = FF(dist, a1, a2);
										}
									}
								en = en * molA.damping_coeff * molB.damping_coeff;
							}
			}
}

results energies_and_forces(const state &molA, const state &molB, double &Lx, double &Ly, double &beta, bool pressure_calc)
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
	double diff_delta = 0.01;


	double U = 0;
	double pressure_X = 0, pressure_Y = 0;
	state molB_clone = molB;

	double dist_x_plus_delta, dist_y_plus_delta;
	double t_U_delta;

	for (int id = -1; id < 2; id++)
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
				if (r2 < min_dist_2)
				{
					en_and_press.energy += E_INF/beta;
					HC_radius = true;
					break;
				}
				if (r2 > min_dist_2 && r2 <= max_dist_2)
				{
					r = sqrt(r2);
					double t_U;
					energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x, dist_y, t_U);
					U += t_U;
					if (pressure_calc)
					{
						t_U_delta = 0;
						molB_clone.x = molB_clone.x - diff_delta;
						dist_x_plus_delta = molB_clone.x - molA.x;
						r2 = dist_x_plus_delta*dist_x_plus_delta + dist_y*dist_y;
						r = sqrt(r2);
						molB_clone.damping_coeff = damping_field(molB_clone.x, Lx);
						molB_clone.ex_field_coeff = external_field(molB_clone.x, Lx);
						energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x_plus_delta, dist_y, t_U_delta);
						pressure_X += -(t_U-t_U_delta)/diff_delta*dist_x;
						molB_clone.x = molB_clone.x + diff_delta;
						molB_clone.damping_coeff = damping_field(molB_clone.x, Lx);
						molB_clone.ex_field_coeff = external_field(molB_clone.x, Lx);

						t_U_delta = 0;
						molB_clone.y = molB_clone.y - diff_delta;
						dist_y_plus_delta = molB_clone.y - molA.y;
						r2 = dist_x*dist_x + dist_y_plus_delta*dist_y_plus_delta;
						r = sqrt(r2);
						energy_calculation(molA, molB_clone, Lx, Ly, beta, r, dist_x, dist_y_plus_delta, t_U_delta);
						pressure_Y += -(t_U-t_U_delta)/diff_delta*dist_y;
						molB_clone.y = molB_clone.y + diff_delta;
					}
				}
		 }
	}
	if (!HC_radius)
	{
		en_and_press.energy = U;
		en_and_press.p_X = pressure_X;
		en_and_press.p_Y = pressure_Y;
		if(en_and_press.energy >= E_INF/beta)
			{
				en_and_press.energy = E_INF/beta;
				en_and_press.p_X = 0;
				en_and_press.p_Y = 0;
			}
	}
	else
  {
    en_and_press.energy = E_INF/beta;
    en_and_press.p_X = 0;
    en_and_press.p_Y = 0;
  }
	return en_and_press;
}

int check_HC(const state &molA, const state &molB, double &Lx, double &Ly)
{
  if ((molA.damping_coeff == 0) || (molB.damping_coeff == 0))
  {
    return 0;
  }
	double r2;	// Distance sqaured
	double dist_x, dist_y;	// Distance between A and B molecules along x and y axies
	double dist_x_2, dist_y_2;


//		double r0 = 7.5877; // Hard core radius in A
	state molB_clone = molB;

	for (int id = -1; id < 2; id++)
	{
     if (HC_radius) {break;}
		 dist_x = 0;
		 dist_y = 0;
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
