double damping_field (double x, double &Lx)
{
	// return 0 - for ideal gas
	// return 1 - for original interactions
	// lambda0 is actually sqrt(lambda0), the same is true for lambdam
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 32.0*abs_x/Lx;
	if(ksi < 5.0){return 1.0;}
		else if(ksi < 7.0){return lambda0 + (1.0 - lambda0)*(7.0 - ksi)*(7.0 - ksi)*(ksi - 4.0)/4.0;}
			else if (ksi < 9.0){return lambda0;}
				else if (ksi < 11.0){return lambdam + (lambda0 - lambdam)*(11.0 - ksi)*(11.0 - ksi)*(ksi - 8.0)/4.0;}
					else {return lambdam;}
	return 1;
}

results external_field (double x, double &Lx)
{
	double abs_x;
	results ex_field_coeff;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 32.0*abs_x/Lx;
	if(ksi < 5.0)
		{
			ex_field_coeff.energy = 0.0;
			ex_field_coeff.p_X = 0.0;
			ex_field_coeff.p_Y = 0.0;
		}
		else if(ksi > 7.0)
					{
						ex_field_coeff.energy = u_m;
						ex_field_coeff.p_X = 0.0;
						ex_field_coeff.p_Y = 0.0;
					}
					else
						{
							ex_field_coeff.energy = u_m*(1.0 - (7.0 - ksi)*(7.0 - ksi)*(ksi - 4.0)/4.0);
							ex_field_coeff.p_X = 12.0*u_m*(7.0 - ksi)*(ksi - 5.0)/Lx;
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
	for (int mol = 0; mol < nPart; mol++)
	{
		double d_lambda = 0.0, d_u_ext = 0;
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}
		double ksi = 32.0*abs_x/Lx; // for damping and external fields from pccp_2022

		// Damping field derivative if the expressions in the damping field definition are in sqrt(lambda0) and for sqrt(lambda)
		// This derivative is correct for damping_field_pccp_2022
		if(ksi < 5.0){d_lambda = 0.0;}
					else if(ksi < 7.0){d_lambda = (1.0 - lambda0)*(ksi - 7.0)*(ksi - 5.0)*(24.0/Lx);}
						else if (ksi < 9.0){d_lambda = 0.0;}
							else if (ksi < 11.0){d_lambda = (lambda0 - lambdam)*(ksi - 11.0)*(ksi - 9.0)*(24.0/Lx);}
								else {d_lambda = 0.0;}

		// External field derivative for external_field_pccp_2022
		if(ksi < 5.0){d_u_ext = 0.0;}
			else if(ksi > 7.0){d_u_ext = 0.0;}
				else{d_u_ext = (24.0/Lx)*u_m*(7.0 - ksi)*(ksi - 5.0);}

				// Calculating the pressure change
				if(coordinates[mol].damping_coeff == 0){EN_AND_PR_counter.p_X += d_u_ext;}
				else
					{
						dp_lambda += - coordinates[mol].en_and_pr.energy*d_lambda/coordinates[mol].damping_coeff;
						dp_ex += d_u_ext;
					}
	}
	EN_AND_PR_counter.p_X = dp_lambda + dp_ex/2.0;
}
