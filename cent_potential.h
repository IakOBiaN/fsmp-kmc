double cent_potential(double x, double Lx)
	{
		double h = 1.0;
		double abs_x = abs(x - Lx/2.0);
		double scale_factor = 1.0;
		if (abs_x/Lx < (Lx/8.0 - scale_factor*max_dist*0.5)/Lx){h = 1.0;}
			else if (abs_x/Lx >= (Lx/8.0 + scale_factor*max_dist*0.5)/Lx){h = 0;}
					else
							{
								double s =  2.0*abs_x/(max_dist*scale_factor) - Lx/(4.0*max_dist*scale_factor);
								h = (2.0 - 3.0*s + s*s*s)/4.0;
							}
		return h;
	}
