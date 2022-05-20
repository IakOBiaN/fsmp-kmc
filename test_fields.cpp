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

int main()
{
	double Lx = 100;
	double abs_x;

	double lambda = 1.0;
	double u_ext = 0;
	double lambda0 = 0.3;
	double u_m = -20;

	cout << "x" << " " << "lambda" << " " << "u_ext" << endl;
	for (double x = 0; x <= Lx; x += 1.0)
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

		cout << x - Lx/2.0 << " " << lambda << " " << u_ext << endl;
	}

	return 0;
}
