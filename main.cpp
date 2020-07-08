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
double energy_QQ;
double p_X, p_Y;
results();      //constructor
results operator+(const results& b) {
         results res;
         res.energy = this->energy + b.energy;
				 res.energy_QQ = this->energy_QQ + b.energy_QQ;
         res.p_X = this->p_X + b.p_X;
         res.p_Y = this->p_Y + b.p_Y;
         return res;
      }
results operator-(const results& b) {
         results res;
         res.energy = this->energy - b.energy;
				 res.energy_QQ = this->energy_QQ - b.energy_QQ;
         res.p_X = this->p_X - b.p_X;
         res.p_Y = this->p_Y - b.p_Y;
         return res;
      }
results operator/(double b) {
         results res;
         res.energy = this->energy/b;
				 res.energy_QQ = this->energy_QQ/b;
         res.p_X = this->p_X/b;
         res.p_Y = this->p_Y/b;
         return res;
      }
results operator*(double b) {
         results res;
         res.energy = this->energy*b;
				 res.energy_QQ = this->energy_QQ*b;
         res.p_X = this->p_X*b;
         res.p_Y = this->p_Y*b;
         return res;
      }
};

//constructor
results::results(void) {
   energy = 0;
	 energy_QQ = 0;
   p_X = 0;
   p_Y = 0;
}

// Class contains the function descrabing the state of the molecule:
// 1) coordinates
// 2) angles
// 3) charge distribution
class state {
public:
double x;
double y;
double q1x_p;
double q1y_p;
double q2x_p;
double q2y_p;
double q3x_p;
double q3y_p;
double q1x_n;
double q1y_n;
double q2x_n;
double q2y_n;
double q3x_n;
double q3y_n;
double phi;
double sin_phi;
double cos_phi;
double energy;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

results EN_AND_PR_counter;											// energy and pressures in the system.
double ACCEPTANCE_RATIO_r[2] = {0, 0};					// 0 - not accepted steps of rotation, 1 - accepted steps of rotation
double ACCEPTANCE_RATIO_m[2] = {0, 0};					// 0 - not accepted steps of move, 1 - accepted steps of move
int BALANCE_STEPS = 100;												// steps for balance statistics
double delta = 1.5;															// MC parameter. Maximal shift of the molecule
double delta_angle = 60.0;    										// MC parameter. Maximal rotation in degrees
double R = 8.3144598;														// Gas constant in J per mol
double N_a = 6.02214e+23;												//	Avogadro constant
const double PI = 3.14159265358979323846;
double density;																	// Actual density of the layer in mkMol per m^2
int N_test;																			// Counter for attempts to insert the test particle in Widom's algorythm
double e_test;																	// Counter for energy change due to the insertion of the test particle

// Minimal (hard core distance) and maximal distance between the molecules
// Hard core radius = 7.5877 A
double min_dist = 7.5877, max_dist = 20.0;

int frame = 0; // For visualization purpose

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "PBC2D.h"
#include "energies_and_forces_exact.h"
//#include "energies_and_forces.h"
#include "initConfigHoneycombTMA.h"
#include "write_xyz_file.h"
#include "PotentialEnergy.h"
#include "Metropolis_iteration.h"
#include "Widom_test.h"

int main()
{
 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 // Set configuration parameters

 double press_X = 0, press_Y = 0, Energy=0, Energy_QQ = 0;
 double persent = 0;

 /////////////////////////////
 // Set the Monte Carlo run //
 /////////////////////////////
 int nPart = 512; // Amount of molecules in the layer
 int nSteps = 500000;            // Total amount of MCS
 int nIter = nSteps * nPart;
 int nStepsEq = 100000;           // MCS for relaxation
 int nIterEq = nStepsEq * nPart;
 double Lx, Ly;  // Linear size of the system in A
 vector <state> coordinates(nPart); // Vector of the molecules coordinates, angles and charges

 ////////////////////////////////////////////////////////////
 //         MC simulation of systems with different N      //
 ////////////////////////////////////////////////////////////

 //Generete an initial distribution of molecules at fixed density
 double temperature = 300.0;
 double deltaT = 10.0;

 initConfigHoneycombTMA(nPart, density, coordinates, Lx, Ly);
 write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, 0, 1, true);

 for(temperature = 300; temperature < 310; temperature += deltaT)
 {

//	vector <double> pressure_avr_X;
//	vector <double> pressure_avr_Y;
	frame = 1;
	write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, frame, 1, false);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////// SYSTEM COUNTERS //////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.energy_QQ = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	press_X = 0;
	press_Y = 0;
	ACCEPTANCE_RATIO_r[0] = 0;  // rejected rotations
	ACCEPTANCE_RATIO_r[1] = 0;  // accepted rotations
	ACCEPTANCE_RATIO_m[0] = 0;  // rejected translations
	ACCEPTANCE_RATIO_m[1] = 0; // accepted translations
	int balanceEq = 0;
	persent = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////

	double beta = 1.0 / (R*temperature);  // Inverse temperature in units of (k_B*T)^-1

/// Write the item for average image
	vector <vector <double> > xy_matrix(6000, vector<double> (6000));
	for(int i = 0; i < 6000; i++){for(int j = 0; j < 6000; j++){xy_matrix[i][j] = 0;}}
// Calculate initial energy
	PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
	cout << "Energy: " << EN_AND_PR_counter.energy/1000.0/nPart << "\t" << "Energy_QQ: " << EN_AND_PR_counter.energy_QQ/1000.0/nPart << "\t" << "P_T: " << (R*temperature*(1.0e+23)*nPart/(Lx*Ly)/N_a)+((EN_AND_PR_counter.p_X + EN_AND_PR_counter.p_Y)/2.0/Lx/Ly*1e23/N_a)<< endl;
	//////////////////////////////////////////////////
	//             Monte Carlo Simulation           //
	//////////////////////////////////////////////////
	double sum_iterations = 0;
	double percent = 0;
	for(int iter = 1; iter <= nIter; iter++)
		{
			percent += 1;
			if(percent > nIter/100.0)
				{
					frame++;
					write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, frame, 1, false);
					cout << int(iter*100.0/nIter) << " %" << endl;
					percent = 0;
				}
			Metropolis_iteration(nPart, Lx, Ly, beta, coordinates);
			if(iter > nIterEq)
				{
					if (iter == nIterEq+1){Energy = 0; Energy_QQ = 0; press_X = 0; press_Y = 0; sum_iterations = 0; N_test = 0; e_test = 0;}
					sum_iterations += 1;
					Energy += EN_AND_PR_counter.energy;
					Energy_QQ += EN_AND_PR_counter.energy_QQ;
					press_X += EN_AND_PR_counter.p_X;
					press_Y += EN_AND_PR_counter.p_Y;
					Widom_test(nPart, coordinates, Lx, Ly, beta, N_test, e_test);
				}
		}

	Energy /= sum_iterations;
	Energy_QQ /= sum_iterations;
	press_X /= sum_iterations;
	press_Y /= sum_iterations;
	press_X *= (1.0/Lx/Ly*1e23)/N_a;
	press_Y *= (1.0/Lx/Ly*1e23)/N_a;
	double mu_ex = log(N_test/(e_test));

	cout << "Energy_MC: " << Energy/1000.0/nPart << " Energe_QQ: " << Energy_QQ/1000.0/nPart << " Pressure: " << (R*temperature*(1.0e+23)*nPart/(Lx*Ly)/N_a) + (press_X + press_Y)/2.0 << " Chem. potential: " << mu_ex << endl;
 }
 return 0;
}
