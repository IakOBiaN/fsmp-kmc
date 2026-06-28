double damping_field (double x, double &Lx)
{
	// return 0 - for ideal gas
	// return 1 - for original interactions
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 32.0*abs_x/Lx - 8.0;
	if(ksi < -3.0){return 1.0;}
		else if(ksi < -1.0){return 1.0 + (1.0 - lambda0)*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;}
			else if (ksi < 1.0){return lambda0;}
				else if (ksi < 3.0){return lambda0*ksi*(ksi*ksi - 6.0*ksi + 9.0)/4.0;}
					else {return 0.0;}
	return 1;
}

results external_field (double x, double &Lx)
{
	double abs_x;
	results ex_field_coeff;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 32.0*abs_x/Lx - 8.0;
	if(ksi < -3.0)
		{
			ex_field_coeff.energy = 0.0;
			ex_field_coeff.p_X = 0.0;
			ex_field_coeff.p_Y = 0.0;
		}
		else if(ksi > -1.0)
					{
						ex_field_coeff.energy = u_m;
						ex_field_coeff.p_X = 0.0;
						ex_field_coeff.p_Y = 0.0;
					}
					else
						{
							ex_field_coeff.energy = -u_m*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;
							ex_field_coeff.p_X = -u_m*24.0*(ksi*ksi + 4.0*ksi + 3.0)/Lx;
							ex_field_coeff.p_Y = 0;
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

	for (int mol = 0; mol < nPart; mol++)
	{
		double d_lambda = 0.0, d_u_ext = 0;
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}
		double ksi = 32.0*abs_x/Lx - 8.0; // for damping and external fields from jpcc_2021

		// Damping field derivative if the expressions in the damping field definition are in sqrt(lambda0) and for sqrt(lambda)
		// This derivative is correct for damping_field_jpcc_2021
		if(ksi < -3.0){d_lambda = 0.0;}
					else if(ksi < -1.0){d_lambda = 6.0*((1.0 - lambda0)/4.0)*((1.0 - lambda0)/4.0)*(ksi*ksi + 4.0*ksi + 3.0)*((ksi*ksi*ksi + 6.0*ksi*ksi + 9.0*ksi) + 4.0/(1 - lambda0))*(32.0/Lx);}
						else if (ksi < 1.0){d_lambda = 0.0;}
							else if (ksi < 3.0){d_lambda = (12.0/Lx)*lambda0*lambda0*ksi*(ksi*ksi - 6.0*ksi + 9.0)*(ksi*ksi - 4.0*ksi + 3.0);}
								else {d_lambda = 0.0;}
			d_lambda = d_lambda/2.0;

		// External field derivative for external_field_jpcc_2021
		if(ksi < -3.0){d_u_ext = 0.0;}
			else if(ksi > -1.0){d_u_ext = 0.0;}
						else{d_u_ext = -(12.0/Lx)*u_m*(ksi*ksi + 4.0*ksi + 3.0);}

		// Calculating the pressure change
		if(coordinates[mol].damping_coeff == 0){EN_AND_PR_counter.p_X += d_u_ext;}
		else
			{
				EN_AND_PR_counter.p_X = EN_AND_PR_counter.p_X - coordinates[mol].en_and_pr.energy*abs(d_lambda)/coordinates[mol].damping_coeff/coordinates[mol].damping_coeff + d_u_ext;
			}
	}
}
