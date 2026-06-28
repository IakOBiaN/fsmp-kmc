double damping_field (double x, double &Lx)
{
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0){return 1.0;}
		else if(ksi > 1.0){return 0;}
			else {return 1.0 - (2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
	return 0;
}

// External field returns a results struct (energy + pressure gradient) to match the
// current program_body, which stores it in state.ex_field_coeff and reads .energy.
results external_field (double x, double &Lx)
{
	double abs_x;
	results ex_field_coeff;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0)
		{
			ex_field_coeff.energy = 0.0;
			ex_field_coeff.p_X = 0.0;
			ex_field_coeff.p_Y = 0.0;
		}
		else if(ksi > 1.0)
			{
				ex_field_coeff.energy = u_m;
				ex_field_coeff.p_X = 0.0;
				ex_field_coeff.p_Y = 0.0;
			}
			else
				{
					ex_field_coeff.energy = u_m*(2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;
					ex_field_coeff.p_X = 12.0*(u_m/Lx)*(1.0 - ksi*ksi);
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

// Extra tangential pressure produced by the gradients of the imposed fields at the
// gas-crystal interface: dp_lambda (damping field) + dp_ex (external field). This is
// the term that closes the mechanical balance  p = p_gas + dp_lambda + dp_ex.
// The derivatives below are d(sqrt(lambda))/dx and d(u_ext)/dx for the jpcc2020
// fields, where both fields use the same ksi = 16*abs_x/Lx - 4 - damping_delta.
void pressure_change_over_interface(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	double dp_lambda = 0, dp_ex = 0;
	for (int mol = 0; mol < nPart; mol++)
	{
		double d_lambda = 0.0, d_u_ext = 0.0;
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}
		double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;

		// d(sqrt(lambda))/dx for the jpcc2020 damping field
		if(ksi < -1.0){d_lambda = 0.0;}
			else if(ksi <= 1.0){d_lambda = (12.0/Lx)*(ksi*ksi - 1.0);}
				else {d_lambda = 0.0;}

		// d(u_ext)/dx for the jpcc2020 external field
		if(ksi < -1.0){d_u_ext = 0.0;}
			else if(ksi <= 1.0){d_u_ext = 12.0*(u_m/Lx)*(1.0 - ksi*ksi);}
				else {d_u_ext = 0.0;}

		// Damping part: molecule interaction energy (minus its external-field energy)
		// times |d(sqrt(lambda))/dx| / sqrt(lambda). External part: direct gradient.
		if(coordinates[mol].damping_coeff == 0){dp_ex += d_u_ext;}
		else
			{
				dp_lambda += - (coordinates[mol].en_and_pr.energy - coordinates[mol].ex_field_coeff.energy)*abs(d_lambda)/coordinates[mol].damping_coeff;
				dp_ex += d_u_ext;
			}
	}
	EN_AND_PR_counter.p_X = dp_lambda + dp_ex/2.0;
}
