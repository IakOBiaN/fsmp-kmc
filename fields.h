double damping_field_jpcc_2020 (double x, double &Lx)
{
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0){return 1.0;}
		else if(ksi > 1.0){return 0;}
			else {return 1.0 - (2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
	return 0;
}

double external_field_jpcc_2020 (double x, double &Lx)
{
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0){return 0.0;}
		else if(ksi > 1.0){return u_m;}
			else {return u_m*(2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
			return 0;
}

//
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
							ex_field_coeff.p_Y = ex_field_coeff.p_X;
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
