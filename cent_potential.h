double cent_potential(double x, double Lx)
	{
		double h;
		double abs_x = abs(x - Lx/2.0);
		if (abs_x/Lx < (Lx/8.0 - max_dist*0.5)/Lx){h = 1.0;}
			else if (abs_x >= (Lx/8.0 + max_dist*0.5)/Lx){h = 0;}
					else
							{
								double s =  2.0*abs_x/max_dist - Lx/(4.0*max_dist);
								h = (2.0 - 3.0*s + s*s*s)/4.0;
							}
		return h;
	}
