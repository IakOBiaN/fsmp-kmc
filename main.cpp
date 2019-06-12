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

bool rosenbluth = false; //kMC NOT WORKING NOW!!!// If rosenbluth = false then Metropolis algorithm works
results EN_AND_PR_counter;                           //energy and pressures in the system.
double ACCEPTANCE_RATIO_r[2] = {0, 0};               //0 - not accepted steps of rotation, 1 - accepted steps of rotation
double ACCEPTANCE_RATIO_m[2] = {0, 0};               //0 - not accepted steps of move, 1 - accepted steps of move
int BALANCE_STEPS = 100;                             //steps for balance statistics
double delta = 1.5;                            //MC parameter. Maximal shift of the molecule
double delta_angle = 60.0;    //MC parameter. Maximal rotation in degrees
double R = 8.3144598;
double N_a = 6.02214e+23;
const double PI  =3.14159265358979323846;
double temperature = 800.0;
double density;                // Actual density of the layer in mkMol per m^2
bool close = false;

// Forcefield for N2-N2 pair
vector <vector <vector <double> > > forcefield;
vector <vector <vector <double> > > energy;
// Minimal (hard core distance) and maximal distance between the molecules
double min_dist,max_dist;
// Delta between neighbor distances in the forcefield in A
double dr;
// Delta between orientation angle of the single molecule
double da;
int frame = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "read_forcefield.h"
#include "PBC2D.h"
#include "energies_and_forces.h"
#include "initConfigHoneycombTMA.h"
#include "initConfigFlowerTMA.h"
#include "initConfig2FlowerTMA.h"
#include "initConfigSuperFlowerTMA.h"
#include "PotentialEnergy.h"
#include "Metropolis_iteration.h"
#include "pressure_balance.h"
#include "layer_map.h"
#include "write_xy_matrix.h"
#include "writeData.h"
#include "write_xyz_file.h"
#include "density_change.h"

int main()
{
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// READ FORCEFIELD FROM FILE ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

  // Fill in the forcefield
  // First dimension is distance
  // Second dimension is angle of first molecule
  // Third dimension is angle of second molecule
  for (int i = 0; i < 551; i++) {
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
  read_forcefield ("Dreiding_ff_TMA_new.dat", forcefield, min_dist,max_dist, dr, da);

//////////////////////////////////////////////////////////////////////////////////////////


 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 // Set configuration parameters

 double Pt = 0;
 double press_X = 0, press_Y = 0, Energy=0;
 double en_2_av = 0;
 double cap_n = 1.0;
 double fluent_capacity = 0;
 double persent = 0;

 /////////////////////////////
 // Set the Monte Carlo run //
 /////////////////////////////
 int nPart = 900; // Amount of molecules in the layer
 int nSteps =700000;            // Total amount of MCS
 int nIter = nSteps * nPart;
 int nStepsEq = 400000;           // MCS for relaxation
 int nIterEq = nStepsEq * nPart;
 double Lx, Ly;  // Linear size of the system in A
 double state_dens = 1.291163; // mkMol of TMA per A^2
 vector <state> coordinates(nPart); // Vector of the molecules coordinates and angles

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Temperature" << "\t" << "Heat.Capacity(reccurent)" << "\t" << "Heat.Capacity" << "\t" << "E_per_molecule" << "\t" << "p_X" << "\t" << "p_Y" << "\t" << "Lx" << "\t" << "Ly" << endl;
 fileOutput.close();

 ////////////////////////////////////////////////////////////
 //         kMC simulation of systems with different N     //
 ////////////////////////////////////////////////////////////

 //for(int nPart = minPart; nPart < maxPart; nPart += stepPart)

 //Generete an initial distribution of molecules at fixed density
initConfigHoneycombTMA(nPart, density, coordinates, Lx, Ly, state_dens);

double deltaT = 100.0;
for(temperature = 400; temperature < 1200; temperature += deltaT)
{
//double delta_rho = 0.1;
//for (density = density; density < 2.6; density += delta_rho)

		// Change the current density of the layer
		//density_change(Lx, Ly, nPart, coordinates);

		 if (temperature > 750){deltaT = 10.0;}
		 write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, 0, 1, true);
     frame = 1;
     write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, frame, 1, false);

     ///////////////////////////////////////////////////////////////////////////////////////////////////
     //////////// SYSTEM COUNTERS //////////////////////////////////////////////////////////////////////
     //////////////////////////////////////////////////////////////////////////////////////////////////
     EN_AND_PR_counter.energy = 0;
     EN_AND_PR_counter.p_X = 0;
     EN_AND_PR_counter.p_Y = 0;
  	 Pt = 0;
     press_X = 0;
     press_Y = 0;
  	 ACCEPTANCE_RATIO_r[0] = 0;  // rejected rotations
  	 ACCEPTANCE_RATIO_r[1] = 0;  // accepted rotations
  	 ACCEPTANCE_RATIO_m[0] = 0;  // rejected translations
  	 ACCEPTANCE_RATIO_m[1] = 0; // accepted translations
  	 double Time = 0; // Total time of the equilibrium run
  	 double Mconf = 0; // Amount of configurations for chemical potential calculation with kMC
  	 double dt = 0;
  	 int balanceEq = 0;
     persent = 0;
     ////////////////////////////////////////////////////////////////////////////////////////////////////

	   double beta = 1.0 / (R*temperature);  // Inverse temperature in units of (k_B*T)^-1

     // Write the item for average image
     vector <vector <double> > xy_matrix(3000, vector<double> (3000));
     for(int i = 0; i < 3000; i++){for(int j = 0; j < 3000; j++){xy_matrix[i][j] = 0;}}

     // Calculate initial energy
		 PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
		 cout << "Energy=" << EN_AND_PR_counter.energy/1000.0/nPart << endl;

     //////////////////////////////////////////////////
     //             Monte Carlo Simulation           //
     //////////////////////////////////////////////////

     for(int iter = 1; iter <= nIter; iter++)
        {

          if((iter%(nPart*500) == 0) || (iter == 1))
          {
           frame++;
           write_xyz_file_TMA (nPart, Lx, Ly, temperature, coordinates, frame, 1, false);
          }
         persent += 1;
         if(persent > nIter/100.0)
         {
           cout << "T = " << temperature << " rho = " << density << " " << int(iter*100.0/nIter) << " %" << endl;persent = 0;
         }

         int trialPart;
         // Choose a particle to be displaced
         // according to the ROSENBLUTH scheme
         // and calculate the duration of the current configuration
         if(!rosenbluth) {Metropolis_iteration(nPart, Lx, Ly, beta, coordinates); dt = 1.0;}     // Make a MC iteration
         //else {trialPart = Rosenbluth_algorithm_simple(nPart, coordinates, dt);}                 // kMC trial particle and dt calculation


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Pressure balance

        balanceEq++;
        if((iter < nIterEq) && (balanceEq > nPart*0.1*BALANCE_STEPS))
        {
            Pt += dt;
            press_X += EN_AND_PR_counter.p_X*dt;
            press_Y += EN_AND_PR_counter.p_Y*dt;
            if(((iter%(BALANCE_STEPS*nPart))==0 && iter != 0) || iter==nIterEq-1)
            {
                press_X /= Pt;
                press_Y /= Pt;
                if (iter < 0.09*nIterEq && iter >= 0.05*nIterEq) { BALANCE_STEPS = 200; }
                if (iter < 0.15*nIterEq && iter >= 0.09*nIterEq) { BALANCE_STEPS = 300; }
                if (iter < 0.25*nIterEq && iter >= 0.15*nIterEq) { BALANCE_STEPS = 500; }
                if (iter < 0.46*nIterEq && iter >= 0.25*nIterEq) { BALANCE_STEPS = 1000; }
                if (iter >= 0.46*nIterEq) { BALANCE_STEPS = 2500; }

                pressure_balance (press_X, press_Y, Lx, Ly, nPart, coordinates, beta);
                Pt = 0;
                press_X = 0;
                press_Y = 0;
                balanceEq = 0;
                ACCEPTANCE_RATIO_r[0] = 0;
                ACCEPTANCE_RATIO_r[1] = 0;
                ACCEPTANCE_RATIO_m[0] = 0;
                ACCEPTANCE_RATIO_m[1] = 0;
            }
        }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

         // Collect the characteristics of interest at equilibrium
         if(iter > nIterEq + nPart * 1000)
           {
            if (iter == nIterEq+1)
            {
                Energy = 0;
                en_2_av = 0;
                press_X = 0;
                press_Y = 0;
                Pt = 0;
            }

            // Update the total time and the amount of configurations for kMC (chemical potential calculation)
            if(rosenbluth) {Time += dt; Mconf++;}

            layer_map_TMA(nPart, coordinates, xy_matrix, Lx, Ly);
            if (cap_n > 1) {fluent_capacity = (cap_n - 1.0)/cap_n*fluent_capacity + (cap_n-1.0)/(cap_n*cap_n)*(Energy/(Pt+1) - EN_AND_PR_counter.energy)*(Energy/(Pt+1) - EN_AND_PR_counter.energy);}

            cap_n++;
            Pt += dt;
            Energy += EN_AND_PR_counter.energy*dt;
            en_2_av += EN_AND_PR_counter.energy*EN_AND_PR_counter.energy*dt;
            press_X += EN_AND_PR_counter.p_X*dt;
            press_Y += EN_AND_PR_counter.p_Y*dt;
           }
         // A new random position is chosen uniformly over the whole volume of the system
         // and update the energies of all molecules in the system
         //if(rosenbluth) {replace_the_trialParticle_and_update_energies(nPart, trialPart, Lx, Ly, beta, coordinates);}
       }

     double mu = 0;
     if(rosenbluth) {mu = log(Mconf/Lx/Ly) - log(Time);}

            press_X /= Pt;
            press_Y /= Pt;

            press_X *= (1.0/Lx/Ly*1e20*1000)/N_a;  //it means p_x_lj = p_x_lj/Lx/Ly/sigma/sigma*1e20*1000 mN/m
            press_Y *= (1.0/Lx/Ly*1e20*1000)/N_a;

            /*press_X = R*temperature*nPart/Ly/Lx + press_X_LJ + press_X_QQ;
            press_Y = R*temperature*nPart/Ly/Lx + press_Y_LJ + press_Y_QQ;*/

            Energy = Energy/Pt;
            en_2_av = en_2_av/Pt;

     // Write the calculated data to a file
     writeData(temperature, fluent_capacity/R/temperature/temperature, (en_2_av-pow(Energy,2))/R/temperature/temperature, Energy/1000.0/nPart, press_X, press_Y, Lx, Ly);

     // Write the xy-matrix
     write_xy_matrix(nPart, Lx, Ly, temperature, xy_matrix);

		 cout << "rho: " << density << " mkMol/m^2 \t" << "mu: " << mu << "\t" << "en: " << Energy/nPart/1000.0 << " kJ/mol" << endl;
     cout << "temp=" << temperature << " P_X=" << press_X << " P_Y=" << press_Y  << endl;

    }
 return 0;
}
