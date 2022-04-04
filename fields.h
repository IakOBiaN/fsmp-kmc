double damping_field (double x, double &Lx)
{
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0){return 1.0;}
		else if(ksi > 1.0){return 0;}
			else {return 1.0 - (2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
}

double external_field (double x, double &Lx)
{
	double abs_x;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 16.0*abs_x/Lx - 4.0 - damping_delta;
	if(ksi < -1.0){return 0.0;}
		else if(ksi > 1.0){return u_m;}
			else {return u_m*(2.0 + 3.0*ksi - ksi*ksi*ksi)/4.0;}
}
