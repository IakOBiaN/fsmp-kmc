#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "random/sfmt.h"
#include "random/sfmt.cpp"
#include <cmath>
#include <float.h>
using namespace std;


double damping_field (double x, double &Lx, double &lambda0)
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


int main()
{
	double Lx = 100;
	double abs_x;
	double dx = -0.1;

	double lambda = 1.0, sqrt_lambda = 1.0, d_log_lambda = 0.0, d_log_lambda_num = 0, d_log_sqrt_lambda = 0, d_log_sqrt_lambda_num = 0;
	double u_ext = 0, d_u_ext = 0.0;
	double lambda0 = 0.3;
	double u_m = -20;

	cout << "x" << " " << "sqrt(lambda)" << " " << "lambda" << " " << "u_ext" << " " << "d_u_ext" << " " << "d_log_lambda" << " " << "d_log_lambda_num" << " " << "d_log_sqrt(lambda)" << " " << "d_log_sqrt(lambda)_num" << endl;
	for (double x = 0; x <= Lx; x += 0.2)
	{
		if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
		double ksi = 32.0*abs_x/Lx - 8.0;

		// Testing the damping field
		// return 0 - for ideal gas
		// return 1 - for original interactions
		if(ksi < -3.0){sqrt_lambda = 1.0;}
			else if(ksi < -1.0){sqrt_lambda = 1.0 + (1.0 - lambda0)*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;}
				else if (ksi < 1.0){sqrt_lambda = lambda0;}
					else if (ksi < 3.0){sqrt_lambda = lambda0*ksi*(ksi*ksi - 6.0*ksi + 9.0)/4.0;}
						else {sqrt_lambda = 0.0;}
		lambda = sqrt_lambda*sqrt_lambda;

		// Testing the external field
		if(ksi < -3.0){u_ext = 0.0;}
			else if(ksi > -1.0){u_ext = u_m;}
						else{u_ext = -u_m*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;}

		// Testing the fields derivatives
		// Damping field derivative
/*		if(ksi < -3.0){d_lambda = 0.0;}
			else if(ksi < -1.0){d_lambda = (24.0/Lx)*(1.0 - lambda0)*(ksi*ksi + 4.0*ksi + 3.0);}
				else if (ksi < 1.0){d_lambda = 0.0;}
					else if (ksi < 3.0){d_lambda = (24.0/Lx)*lambda0*(ksi*ksi - 4.0*ksi + 3.0);}
						else {d_lambda = 0.0;}
*/
		// Derivative of the sqrt(damping field)
		if(ksi < -3.0){d_log_sqrt_lambda = 0.0;}
			else if(ksi < -1.0){d_log_sqrt_lambda = (96.0/Lx)*((ksi*ksi + 4*ksi + 3)/(ksi*ksi*ksi + 6*ksi*ksi + 9*ksi + (4/(1-lambda0))));}
				else if (ksi < 1.0){d_log_sqrt_lambda = 0.0;}
					else if (ksi < 3.0){d_log_sqrt_lambda = (96.0/Lx)*((ksi*ksi - 4.0*ksi + 3.0)/(ksi*ksi*ksi - 6*ksi*ksi +9*ksi));}
						else {d_log_sqrt_lambda = 0.0;}

		// Derivative of the damping field
		if(ksi < -3.0){d_log_lambda = 0.0;}
			else if(ksi < -1.0){d_log_lambda = (48.0/Lx)*(2.0*ksi*(ksi + 3.0) + (ksi*ksi + 6.0*ksi + 9.0))/(ksi*(ksi*ksi + 6.0*ksi + 9.0) + (4.0/(1-lambda0)));}
				else if (ksi < 1.0){d_log_lambda = 0.0;}
					else if (ksi < 3.0){d_log_lambda = (64.0/Lx)*(2.0*ksi*(ksi - 3.0) + (ksi*ksi - 6.0*ksi + 9.0))/(ksi*(ksi*ksi - 6.0*ksi + 9.0));}
						else {d_log_lambda = 0.0;}

		// Numeric derivatives of lambda
		double k;
		double dop_sqrt_lambda_plus_dx = 0, dop_sqrt_lambda_min_dx = 0;
		if (x < Lx/2.0){k = -1;} else {k = 1;}
		if(ksi < 3)
			{
				dop_sqrt_lambda_plus_dx = damping_field(x+dx/2.0, Lx, lambda0);
				dop_sqrt_lambda_min_dx = damping_field(x-dx/2.0, Lx, lambda0);
				d_log_sqrt_lambda_num = (log(dop_sqrt_lambda_plus_dx) - log(dop_sqrt_lambda_min_dx))/dx*k;
				d_log_lambda_num = (log(dop_sqrt_lambda_plus_dx*dop_sqrt_lambda_plus_dx) - log(dop_sqrt_lambda_min_dx*dop_sqrt_lambda_min_dx))/dx*k;
			}
			else {d_log_sqrt_lambda_num = 0; d_log_lambda_num = 0;}

		// External field derivative
		if(ksi < -3.0){d_u_ext = 0.0;}
			else if(ksi > -1.0){d_u_ext = 0.0;}
						else{d_u_ext = -(24.0/Lx)*u_m*(ksi*ksi + 4.0*ksi + 3.0);}

		cout << x - Lx/2.0 << " " << sqrt_lambda << " " << lambda << " " << u_ext << " " << d_u_ext << " " << d_log_lambda << " " << d_log_lambda_num << " " << d_log_sqrt_lambda << " " << d_log_sqrt_lambda_num << endl;
	}

	return 0;
}
