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
double phi;
double sin_phi;
double cos_phi;
double damping_coeff;
double ex_field_coeff;
double stat_weight;
double energy;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

results EN_AND_PR_counter;											// energy and pressures in the system.
double ACCEPTANCE_RATIO_r[2] = {0, 0};					// 0 - not accepted steps of rotation, 1 - accepted steps of rotation
double ACCEPTANCE_RATIO_m[2] = {0, 0};					// 0 - not accepted steps of move, 1 - accepted steps of move
int BALANCE_STEPS = 100;												// steps for balance statistics
double delta = 2.5;															// MC parameter. Maximal shift of the molecule
double delta_angle = 60.0;    										// MC parameter. Maximal rotation in degrees
double R = 8.31446261815324;														// Gas constant in J per mol
double N_a = 6.02214076e+23;												//	Avogadro constant
const double PI = 3.14159265358979323846;
double density;																	// Actual density of the layer in mkMol per m^2
double state_dens, state_Ly;
int N_test;																			// Counter for attempts to insert the test particle in Widom's algorythm
double e_test;																	// Counter for energy change due to the insertion of the test particle
double damping_delta = 0;												// Small parameter that elongates the damping field along the Lx dimension
double lambda0 = 0.3;
double u_m = -500.0;													// Parameter of the external field
bool HC_radius = false;                         // Is we inside hard core radius (min_dist)?

// Minimal (hard core distance) and maximal distance between the molecules
// Hard core radius = 7.5877 A

// Forcefield for N2-N2 pair
vector <vector <vector <double> > > forcefield;
vector <vector <vector <double> > > energy;
// Minimal (hard core distance) and maximal distance between the molecules
double min_dist,max_dist;
// Delta between neighbor distances in the forcefield in A
double dr;
// Delta between orientation angle of the single molecule
double da;

//double min_dist = 7.5877, max_dist = 11.052*5.0;

int frame = 0; // For visualization purpose

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "read_forcefield.h"
#include "PBC2D.h"
#include "fields.h"
#include "energies_and_forces_numerical_Dreiding_TMA.h"
//#include "energies_and_forces_numerical_simple_model.h"
#include "initConfigHoneycombTMA_elongated_cell.h"
//#include "initConfigHoneycombTMA.h"
//#include "initConfigFlowerTMA.h"
//#include "initConfigSuperFlowerTMA.h"
#include "write_xyz_file.h"
#include "PotentialEnergy.h"
#include "Metropolis_iteration.h"
#include "pressure_balance.h"
//#include "Widom_test.h"
#include "block_error.h"
#include "density_to_Ly.h"

int main()
{
 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 //////////////////////////////////////////////////////////////////////////////////////
 //////////////////////////// READ FORCEFIELD FROM FILE ///////////////////////////////
 /////////////////////////////////////////////////////////////////////////////////////

	// Fill in the forcefield
	// First dimension is distance
	// Second dimension is angle of first molecule
	// Third dimension is angle of second molecule
	for (int i = 0; i < 2385; i++) {
		vector< vector<double> > mat; // Create an empty matrix
			for (int j = 0; j < 361; j++) {
				vector<double> row; // Create an empty row
					for (int k =0; k <361; k++) {
						row.push_back(0);
					}
					mat.push_back(row); // Add an element (column) to the row
			}
			forcefield.push_back(mat); // Add the row to the main vector
	}
   // Read the forcefield from "forcefield.dat"
   cout << "read forcefield.dat file" << endl;
   read_forcefield ("simplified_model_num_potential_r_7.58_5524_002_phi_1.dat", forcefield, min_dist, max_dist, dr, da);

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Temperature, K" << "\t" << "Density, mkmol/m2" << "\t" << "Lx, A" << "\t" << "Ly, A" << "\t" << "Energy, kJ/mol" << "\t" << "Energy SD" << "\t" << "Pressure, mN/m" << "\t" << "Pressure SD" << "\t" << "P_ex" << "\t" << "P_x" << "\t" << "P_y" << endl;
 fileOutput.close();

 // Set configuration parameters

 double press_X = 0, press_Y = 0, Energy=0, Energy_QQ = 0;
 double persent = 0, AR_r, AR_m;

 /////////////////////////////
 // Set the Monte Carlo run //
 /////////////////////////////
 int nPart = 640; // Honeycomb
// int nPart = 864; // Flower-1
// int nPart = 450; // Superflower
 int nSteps = 15000;            // Total amount of MCS
 int nIter = nSteps * nPart;
 int nStepsEq = 10000;           // MCS for relaxation
 int nIterEq = nStepsEq * nPart;
 double Lx, Ly;  // Linear size of the system in A
 vector <state> coordinates(nPart*4); // Vector of the molecules coordinates, angles and charges
 vector <double> pressure_stat(nIter - nIterEq);
 vector <double> energy_stat(nIter - nIterEq);
 vector <double> energy_QQ_stat(nIter - nIterEq);

 ////////////////////////////////////////////////////////////
 //         MC simulation of systems with different N      //
 ////////////////////////////////////////////////////////////

 //Generete an initial distribution of molecules at fixed density
 double temperature = 400.0;
 double deltaT = 2000.0;
 state_dens = 1.283; // Density that you want in mkmol/m2

// for(state_dens = 1.50; state_dens < 1.60; state_dens += 0.005)
for(temperature = 400; temperature <= 2000; temperature += deltaT)
 {
	initConfigHoneycombTMA_elongated_cell(nPart, density, coordinates, Lx, Ly, state_dens);
//	initConfigHoneycombTMA(nPart, density, coordinates, Lx, Ly, density_to_Ly(nPart, state_dens));
//	initConfigFlowerTMA(nPart, density, coordinates, Lx, Ly, density_to_Ly(nPart, state_dens));
//	initConfigSuperFlowerTMA(nPart, density, coordinates, Lx, Ly, density_to_Ly_SF(nPart, state_dens));
	write_xyz_file_TMA (nPart, density, Lx, Ly, temperature, coordinates, 0, 1, true);
	frame = 1;
	write_xyz_file_TMA (nPart, density, Lx, Ly, temperature, coordinates, frame, 1, false);

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
	cout << "Density: " << density << "\t" << " Energy: " << EN_AND_PR_counter.energy/1000.0/nPart << "\t" << " P: " << (R*temperature*(1.0e+23)*nPart/(Lx*Ly)/N_a)+((EN_AND_PR_counter.p_X + EN_AND_PR_counter.p_Y)/2.0/Lx/Ly*1e23/N_a)<< endl;
	cout << "P_X: " << (EN_AND_PR_counter.p_X/Lx/Ly*1e23/N_a) << "\t" << "P_Y: " << (EN_AND_PR_counter.p_Y/Lx/Ly*1e23/N_a) <<  endl;

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
					write_xyz_file_TMA (nPart, density, Lx, Ly, temperature, coordinates, frame, 1, false);
					cout << int(iter*100.0/nIter) << " %" << endl;
					percent = 0;
				}
			Metropolis_iteration(nPart, Lx, Ly, beta, coordinates);
/*
				balanceEq++;
				BALANCE_STEPS = 100;
				if((iter < nIterEq) && (balanceEq > nPart*0.1*BALANCE_STEPS))
					{
						sum_iterations += 1;
						press_X += EN_AND_PR_counter.p_X;
						press_Y += EN_AND_PR_counter.p_Y;

						if(((iter%(BALANCE_STEPS*nPart))==0 && iter != 0) || iter==nIterEq-1)
							{
								press_X /= sum_iterations;
								press_Y /= sum_iterations;
								if (iter < 0.09*nIterEq && iter >= 0.05*nIterEq) { BALANCE_STEPS = 200; }
								if (iter < 0.15*nIterEq && iter >= 0.09*nIterEq) { BALANCE_STEPS = 300; }
								if (iter < 0.25*nIterEq && iter >= 0.15*nIterEq) { BALANCE_STEPS = 500; }
								if (iter < 0.46*nIterEq && iter >= 0.25*nIterEq) { BALANCE_STEPS = 1000; }
								if (iter >= 0.46*nIterEq) { BALANCE_STEPS = 2500; }
								zero_pressure_balance(press_X, press_Y, Lx, Ly, nPart, coordinates, beta);
//								pressure_balance(press_X, press_Y, Lx, Ly, nPart, coordinates, beta);

								AR_r = ACCEPTANCE_RATIO_r[1]/(ACCEPTANCE_RATIO_r[0]+ACCEPTANCE_RATIO_r[1]);
								AR_m = ACCEPTANCE_RATIO_m[1]/(ACCEPTANCE_RATIO_m[0]+ACCEPTANCE_RATIO_m[1]);
								cout << "AR_m: " << AR_m << " delta: " << delta << " AR_r: " << AR_r << " delta_ang: " << delta_angle << endl;
								if (AR_r < 0.25 && delta_angle > 5.0)
								{delta_angle -= 1.0;}
								if (AR_r > 0.3 && delta_angle < 120.0)
								{delta_angle += 1.0;}
								if (AR_m < 0.25 && delta > 0.05)
								{delta -= 0.02;}
								if (AR_m > 0.3 && delta < 0.8)
								{delta += 0.02;}

								sum_iterations = 0;
								press_X = 0;
								press_Y = 0;
								balanceEq = 0;
								ACCEPTANCE_RATIO_r[0] = 0;
								ACCEPTANCE_RATIO_r[1] = 0;
								ACCEPTANCE_RATIO_m[0] = 0;
								ACCEPTANCE_RATIO_m[1] = 0;
							}
					}
*/
			if(iter > nIterEq)
				{
					if (iter == nIterEq+1){Energy = 0; press_X = 0; press_Y = 0; sum_iterations = 0; N_test = 0; e_test = 0;}
					sum_iterations += 1;
					Energy += EN_AND_PR_counter.energy;
					press_X += EN_AND_PR_counter.p_X;
					press_Y += EN_AND_PR_counter.p_Y;
					energy_stat[sum_iterations] = EN_AND_PR_counter.energy;
					pressure_stat[sum_iterations] = (EN_AND_PR_counter.p_X + EN_AND_PR_counter.p_Y)/2.0;
//					Widom_test(nPart, coordinates, Lx, Ly, beta, N_test, e_test);
				}
		}

	Energy /= sum_iterations;
	press_X /= sum_iterations;
	press_Y /= sum_iterations;
	press_X *= (1.0/Lx/Ly*1e23)/N_a;
	press_Y *= (1.0/Lx/Ly*1e23)/N_a;
//	double mu_ex = log(N_test/(e_test));

/////////// Block Error Calculation ////////////
	double energy_error = block_error_calculation(energy_stat, sum_iterations)/1000.0/nPart;
	double pressure_error = block_error_calculation(pressure_stat, sum_iterations)*(1.0/Lx/Ly*1e23)/N_a;


	cout << "Density: " << density << " Lx: " << Lx << " Ly: " << Ly << endl;
	cout << "T: " << temperature << " Energy_MC: " << Energy/1000.0/nPart << " energy_error: " << energy_error << endl;
	cout << "Pressure: " << (R*temperature*(1.0e+23)*nPart/(Lx*Ly)/N_a) + (press_X + press_Y)/2.0 << " pressure_error: " << pressure_error << endl;
	cout << "P_ex_MC: " << (press_X + press_Y)/2.0 << " P_ex_MC_X: " << press_X << " P_ex_MC_Y: " << press_Y << endl;

	ofstream fileOutput("statistics.dat", ios_base::app);
	fileOutput  << temperature << "\t" << density << "\t" << Lx << "\t" << Ly << "\t" << Energy/1000.0/nPart << "\t" << energy_error << "\t" << (R*temperature*(1.0e+23)*nPart/(Lx*Ly)/N_a) + (press_X + press_Y)/2.0 << "\t" << pressure_error << "\t" << (press_X + press_Y)/2.0 << "\t" << press_X << "\t" << press_Y << endl;
	fileOutput.close();

 }
 return 0;
}
