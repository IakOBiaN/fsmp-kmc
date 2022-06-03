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
	double d_lambda_num = 0, d_log_lambda_num = 0;

	double lambda = 1.0, d_lambda = 0.0, d_log_lambda = 0;
	double u_ext = 0, d_u_ext = 0.0;
	double lambda0 = 0.3;
	double u_m = -20;

	cout << "x" << " " << "lambda" << " " << "u_ext" << " " << "d_lambda" << " " << "d_log_lambda" << " " << "d_u_ext" << " " << "d_lambda_num" << " " << "d_log_lambda_num" << endl;
	for (double x = 0; x <= Lx; x += 0.2)
	{
		if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
		double ksi = 32.0*abs_x/Lx - 8.0;

		// Testing the damping field
		// return 0 - for ideal gas
		// return 1 - for original interactions
		if(ksi < -3.0){lambda = 1.0;}
			else if(ksi < -1.0){lambda = 1.0 + (1.0 - lambda0)*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;}
				else if (ksi < 1.0){lambda = lambda0;}
					else if (ksi < 3.0){lambda = lambda0*ksi*(ksi*ksi - 6.0*ksi + 9.0)/4.0;}
						else {lambda = 0.0;}

		// Testing the external field
		if(ksi < -3.0){u_ext = 0.0;}
			else if(ksi > -1.0){u_ext = u_m;}
						else{u_ext = -u_m*ksi*(ksi*ksi + 6.0*ksi + 9.0)/4.0;}


		// Testing the fields derivatives
		// Damping field derivative
		if(ksi < -3.0){d_lambda = 0.0;}
			else if(ksi < -1.0){d_lambda = (24.0/Lx)*(1.0 - lambda0)*(ksi*ksi + 4.0*ksi + 3.0);}
				else if (ksi < 1.0){d_lambda = 0.0;}
					else if (ksi < 3.0){d_lambda = (24.0/Lx)*lambda0*(ksi*ksi - 4.0*ksi + 3.0);}
						else {d_lambda = 0.0;}

		// Log damping field derivative
		if(ksi < -3.0){d_log_lambda = 0.0;}
			else if(ksi < -1.0){d_log_lambda = (96.0/Lx)*((ksi*ksi + 4*ksi + 3)/(ksi*ksi*ksi + 6*ksi*ksi +9*ksi +(4/(1-lambda0))));}
				else if (ksi < 1.0){d_log_lambda = 0.0;}
					else if (ksi < 3.0){d_log_lambda = (96.0/Lx)*((ksi*ksi - 4.0*ksi + 3.0)/(ksi*ksi*ksi - 6*ksi*ksi +9*ksi));}
						else {d_log_lambda = 0.0;}

		// Numeric derivatives of lambda
		double k;
		if (x < Lx/2.0){k = -1;} else {k = 1;}
		d_lambda_num = (damping_field(x+dx/2.0, Lx, lambda0) - damping_field(x-dx/2.0, Lx, lambda0))/dx*k;
		if(ksi < 3){d_log_lambda_num = (log(damping_field(x+dx/2.0, Lx, lambda0)) - log(damping_field(x-dx/2.0, Lx, lambda0)))/dx*k;}
			else {d_log_lambda_num = 0;}

		// External field derivative
		if(ksi < -3.0){d_u_ext = 0.0;}
			else if(ksi > -1.0){d_u_ext = 0.0;}
						else{d_u_ext = -(24.0/Lx)*u_m*(ksi*ksi + 4.0*ksi + 3.0);}

		cout << x - Lx/2.0 << " " << lambda << " " << u_ext << " " << d_lambda << " " << d_log_lambda << " " << d_u_ext << " " << d_lambda_num << " " << d_log_lambda_num << endl;
	}

	return 0;
}
