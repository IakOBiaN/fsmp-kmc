void pressure_change_over_interface(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	double d_lambda = 0.0, d_log_lambda = 0, d_u_ext = 0;

	for (int mol = 0; mol < nPart; mol++)
	{
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}
		double ksi = 32.0*abs_x/Lx - 8.0;

				// Derivative of the sqrt(damping field)
/*				if(ksi < -3.0){d_log_lambda = 0.0;}
					else if(ksi < -1.0){d_log_lambda = (96.0/Lx)*((ksi*ksi + 4*ksi + 3)/(ksi*ksi*ksi + 6*ksi*ksi +9*ksi +(4/(1-lambda0))));}
						else if (ksi < 1.0){d_log_lambda = 0.0;}
							else if (ksi < 3.0){d_log_lambda = (96.0/Lx)*((ksi*ksi - 4.0*ksi + 3.0)/(ksi*ksi*ksi - 6*ksi*ksi +9*ksi));}
								else {d_log_lambda = 0.0;}
*/
				// Derivative of the damping field
				if(ksi < -3.0){d_log_lambda = 0.0;}
					else if(ksi < -1.0){d_log_lambda = (48.0/Lx)*(2.0*ksi*(ksi + 3.0) + (ksi*ksi + 6.0*ksi + 9.0))/(ksi*(ksi*ksi + 6.0*ksi + 9.0) + (4.0/(1-lambda0)));}
						else if (ksi < 1.0){d_log_lambda = 0.0;}
							else if (ksi < 3.0){d_log_lambda = (64.0/Lx)*(2.0*ksi*(ksi - 3.0) + (ksi*ksi - 6.0*ksi + 9.0))/(ksi*(ksi*ksi - 6.0*ksi + 9.0));}
								else {d_log_lambda = 0.0;}

				// External field derivative
				if(ksi < -3.0){d_u_ext = 0.0;}
					else if(ksi > -1.0){d_u_ext = 0.0;}
								else{d_u_ext = -(24.0/Lx)*u_m*(ksi*ksi + 4.0*ksi + 3.0);}

				// Calculating the pressure change
				EN_AND_PR_counter.p_X +=  -coordinates[mol].en_and_pr.energy*d_log_lambda + d_u_ext;
//				cout << "du_ext: " << d_u_ext << " u: " << coordinates[mol].en_and_pr.energy << " d_log_lambda: " << d_log_lambda << " " << "dp: " << EN_AND_PR_counter.p_X << endl;
	}
}
