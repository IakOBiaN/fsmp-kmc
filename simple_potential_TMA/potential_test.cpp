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


// Random number generator
int seed = (int)time(0);
CRandomSFMT0 RanGen(seed);

class results {
public:
double energy;
double p_X, p_Y;
results();      //constructor
results operator+(const results& b) {
         results res;
         res.energy = this->energy + b.energy;
         res.p_X = this->p_X + b.p_X;
         res.p_Y = this->p_Y + b.p_Y;
         return res;
      }
results operator-(const results& b) {
         results res;
         res.energy = this->energy - b.energy;
         res.p_X = this->p_X - b.p_X;
         res.p_Y = this->p_Y - b.p_Y;
         return res;
      }
results operator/(double b) {
         results res;
         res.energy = this->energy/b;
         res.p_X = this->p_X/b;
         res.p_Y = this->p_Y/b;
         return res;
      }
results operator*(double b) {
         results res;
         res.energy = this->energy*b;
         res.p_X = this->p_X*b;
         res.p_Y = this->p_Y*b;
         return res;
      }
};

//constructor
results::results(void) {
   energy = 0;
   p_X = 0;
   p_Y = 0;
}

// Class contains the function descrabing the state of the molecule:
// 1) coordinates
// 2) angles
class state {
public:
double x;
double y;
double phi;
double sin_phi;
double cos_phi;
double energy;
double mob;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

results en_and_press;
double R = 8.3144598;
double N_a = 6.02214e+23;
const double PI = 3.14159265358979323846;
double temperature = 300.0;
double beta = 1.0 / (R*temperature);
int nPart = 2;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "PBC2D.h"
#include "coarse_grain_potential.h"

int main()
{
	double Lx = 300.0;  // Linear size of the system in A
	double Ly = 300.0;  // Linear size of the system in A
	state molA;
	state molB;
	molA.x = PBC2D(Lx, 10.0);
	molA.y = PBC2D(Ly, 10.0);
	molA.phi = 0.0;
	molA.sin_phi = sin(molA.phi/180.0*PI);
	molA.cos_phi = cos(molA.phi/180.0*PI);

	stringstream name;
	name <<  "potential.dat";
	ofstream fileOutput(name.str().c_str(), ios_base::trunc);
	fileOutput << "Distance" << "\t" << "molA angle" << "\t" << "molB angle" << "\t" << "u" << endl;

	for(double x = 20.0; x <= 60.02; x += 0.02)
//	for(double phi = 0.0; phi <= 60.0; phi += 1.0)
	{
//		molB.x = PBC2D(Lx, 21.48);
		molB.x = PBC2D(Lx, x);
		molB.y = PBC2D(Ly, 10.0);
		molB.phi = 60.0;
//		molB.phi = phi;
		molB.sin_phi = sin(molB.phi/180.0*PI);
		molB.cos_phi = cos(molB.phi/180.0*PI);
		double r = sqrt((molB.x - molA.x)*(molB.x - molA.x) + (molB.y - molA.y)*(molB.y - molA.y));
		en_and_press = coarse_grain_potential(molA, molB, Lx, Ly, false);
		fileOutput << r << "\t" << molA.phi << "\t" << molB.phi << "\t" << en_and_press.energy << endl;
	}
	cout << " Done! " << endl;
	fileOutput.close();
}
