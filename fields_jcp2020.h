double damping_field (double x, double &Lx)
{
	// return 0 - for ideal gas
	// return 1 - for original interactions
	// lambda0 is actually sqrt(lambda0), the same is true for lambdam

	double abs_x;
	double sqrt_lambda = 0;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - delta_damp;
	if(ksi < -1){sqrt_lambda = 1.0;}
		else if(ksi <= 1){sqrt_lambda = 1.0 - (2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
			else {sqrt_lambda = 0;} // ksi > 1
	return sqrt_lambda;
}

results external_field (double x, double &Lx)
{
	double abs_x;
	results ex_field_coeff;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0;
	if(ksi < -1.0)
		{
			ex_field_coeff.energy = 0.0;
			ex_field_coeff.p_X = 0.0;
			ex_field_coeff.p_Y = 0.0;
		}
		else if(ksi <= 1.0)
					{
						ex_field_coeff.energy = u_m*(2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;
						ex_field_coeff.p_X = 12.0*(u_m/Lx)*(1.0 - ksi*ksi);
						ex_field_coeff.p_Y = 0.0;
					}
					else
						{
							ex_field_coeff.energy = u_m;
							ex_field_coeff.p_X = 0.0;
							ex_field_coeff.p_Y = 0.0;
						}
	return ex_field_coeff;
}

double weights_for_central_cell (double x, double &Lx)
{
  double abs_x;
  if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
  double s = 32.0*abs_x/Lx - 4.0;
  if (s <= -1){return 1.0;}
    else if(s < 1.0){return (2.0 - 3.0*s + s*s*s)/4.0;}
  return 0;
}

void pressure_change_over_interface(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	double dp_lambda = 0, dp_ex = 0;
	double ksi = 0;
	for (int mol = 0; mol < nPart; mol++)
	{
		double d_lambda = 0.0, d_u_ext = 0;
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}

		// Damping field derivative if the expressions in the damping field definition are in sqrt(lambda0)
		// This derivative is correct for damping_field_jcp_2020
		ksi = 16.0*abs_x/Lx - 4.0 - delta_damp;
		if(ksi < -1){d_lambda = 0;}
			else if(ksi <= 1){d_lambda = (12.0/Lx)*(ksi*ksi - 1);}
				else {d_lambda = 0;} // ksi > 1

		// External field derivative for external_field_pccp_2022
		ksi = 16.0*abs_x/Lx - 4.0;
		if(ksi < -1.0){d_u_ext = 0.0;}
			else if(ksi <= 1.0){d_u_ext = 12.0*(u_m/Lx)*(1.0 - ksi*ksi);}
				else{d_u_ext = 0;}

				// Calculating the pressure change
				if(coordinates[mol].damping_coeff == 0){dp_ex += d_u_ext;}
				else
					{
						dp_lambda += - coordinates[mol].en_and_pr.energy*abs(d_lambda)/coordinates[mol].damping_coeff;
						dp_ex += d_u_ext;
					}
	}
//	cout << "dp_lambda: " << dp_lambda*1e23/Ly/N_a << " dp_u: " << dp_ex*1e23/Ly/N_a/2.0 << endl;
	EN_AND_PR_counter.p_X = dp_lambda + dp_ex/2.0;
}
