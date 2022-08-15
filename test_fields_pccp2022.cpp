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

// Fields from PCCP 2022
double damping_field (double x, double &Lx, double &lambda0, double &lambdam)
{
	// This function returned sqrt(lambda) as it defined in PCCP 2022
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

double external_field (double x, double &Lx, double &u_m)
{
	double abs_x;
	double ex_field_coeff;
	if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}
	double ksi = 32.0*abs_x/Lx;
	if(ksi < 5.0){ex_field_coeff = 0.0;}
		else if(ksi > 7.0){ex_field_coeff = u_m;}
					else{ex_field_coeff = u_m*(1.0 - (7.0 - ksi)*(7.0 - ksi)*(ksi - 4.0)/4.0);}
	return ex_field_coeff;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
	double Lx = 500;
	double abs_x;
	double dx = -0.0001;

	double lambda = 1.0, sqrt_lambda = 1.0, d_lambda_an = 0, d_lambda_num = 0;
	double u_ext = 0, d_u_ext = 0.0, d_u_ext_num = 0.0;
	double lambda0 = 0.3;
	double lambdam = 0.0;
	double u_m = -20;

	cout << "x" << " " << "u_ext" << " " << "d_u_ext_an" << " " << "d_u_ext_num" << " " << "sqrt_lambda" << " " << "d_lambda_an" << " " << "d_lambda_num" << endl;

	for (double x = 0; x <= Lx; x += 1.0)
	{
		if (x > Lx/2.0){abs_x = x - Lx/2.0;} else {abs_x = Lx/2.0 - x;}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Fields from PCCP 2022

		double ksi = 32.0*abs_x/Lx;

		// sqrt(lambda) as it defined in pccp_2022
		if(ksi < 5.0){sqrt_lambda = 1.0;}
			else if(ksi < 7.0){sqrt_lambda = lambda0 + (1.0 - lambda0)*(7.0 - ksi)*(7.0 - ksi)*(ksi - 4.0)/4.0;}
				else if (ksi < 9.0){sqrt_lambda = lambda0;}
					else if (ksi < 11.0){sqrt_lambda = lambdam + (lambda0 - lambdam)*(11.0 - ksi)*(11.0 - ksi)*(ksi - 8.0)/4.0;}
						else {sqrt_lambda = lambdam;}

		// External field pccp_2022
		if(ksi < 5.0){u_ext = 0.0;}
			else if(ksi > 7.0){u_ext = u_m;}
						else{u_ext = u_m*(1.0 - (7.0 - ksi)*(7.0 - ksi)*(ksi - 4.0)/4.0);}

		// Testing the fields derivatives

		// This derivative d(sqrt(lambda))/dx for damping_field_pccp_2022
/*		if(ksi < 5.0){d_lambda_an = 0.0;}
					else if(ksi < 7.0){d_lambda_an = (7.0 - ksi)*(ksi - 5.0)*(24.0/Lx);}
						else if (ksi < 9.0){d_lambda_an = 0.0;}
							else if (ksi < 11.0){d_lambda_an = (11.0 - ksi)*(ksi - 9.0)*(24.0/Lx);}
								else {d_lambda_an = 0.0;}
*/
		// CORRECTED derivative d(sqrt(lambda))/dx for damping_field_pccp_2022
		if(ksi < 5.0){d_lambda_an = 0.0;}
					else if(ksi < 7.0){d_lambda_an = (1.0 - lambda0)*(ksi - 7.0)*(ksi - 5.0)*(24.0/Lx);}
						else if (ksi < 9.0){d_lambda_an = 0.0;}
							else if (ksi < 11.0){d_lambda_an = (lambda0 - lambdam)*(ksi - 11.0)*(ksi - 9.0)*(24.0/Lx);}
								else {d_lambda_an = 0.0;}

		// External field derivative pccp_2022
		if(ksi < 5.0){d_u_ext = 0.0;}
			else if(ksi > 7.0){d_u_ext = 0.0;}
				else{d_u_ext = (24.0/Lx)*u_m*(7.0 - ksi)*(ksi - 5.0);}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Numeric derivatives of lambda
		double k;
		double lambda_plus_dx = 0, lambda_min_dx = 0, u_ext_plus_dx = 0, u_ext_min_dx = 0;
		if (x < Lx/2.0){k = -1;} else {k = 1;}

		lambda_plus_dx = damping_field(x+dx/2.0, Lx, lambda0, lambdam);
		lambda_min_dx = damping_field(x-dx/2.0, Lx, lambda0, lambdam);

		// dlambda/dx if damping field is defined as sqrt(lambda)
		d_lambda_num = (lambda_plus_dx - lambda_min_dx)/dx*k;

		u_ext_plus_dx = external_field(x+dx/2.0, Lx, u_m);
		u_ext_min_dx = external_field(x-dx/2.0, Lx, u_m);
		d_u_ext_num = (u_ext_plus_dx - u_ext_min_dx)/dx*k;

		cout << x - Lx/2.0 << " " << u_ext << " " << d_u_ext << " " << d_u_ext_num << " " << sqrt_lambda << " " << d_lambda_an << " " << d_lambda_num << endl;
	}

	return 0;
}
