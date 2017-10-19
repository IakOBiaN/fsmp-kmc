#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <time.h>
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

class pressure {
public:
double X_LJ;
double X_QQ;
double Y_LJ;
double Y_QQ;
pressure operator+(pressure& b) {
         pressure press;
         press.X_LJ = this->X_LJ + b.X_LJ;
         press.X_QQ = this->X_QQ + b.X_QQ;
         press.Y_LJ = this->Y_LJ + b.Y_LJ;
         press.Y_QQ = this->Y_QQ + b.Y_QQ;
         return press;
      }
pressure operator-(pressure& b) {
         pressure press;
         press.X_LJ = this->X_LJ - b.X_LJ;
         press.X_QQ = this->X_QQ - b.X_QQ;
         press.Y_LJ = this->Y_LJ - b.Y_LJ;
         press.Y_QQ = this->Y_QQ - b.Y_QQ;
         return press;
      }
};

class results {
public:
double energy;
pressure p;
results();      //constructor
results operator+(results& b) {
         results res;
         res.energy = this->energy + b.energy;
         res.p = this->p + b.p;
         return res;
      }
results operator-(results& b) {
         results res;
         res.energy = this->energy - b.energy;
         res.p = this->p - b.p;
         return res;
      }
};

//constructor
results::results(void) {
   energy = 0;
   p.X_LJ = 0;
   p.X_QQ = 0;
   p.Y_LJ = 0;
   p.Y_QQ = 0;
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

bool rosenbluth = false; //kMC NOT WORKING NOW!!!// If rosenbluth = false then Metropolis algorithm works
bool energy_QQ_exact = false;
bool pressure_QQ_exact = false;
results EN_AND_PR_counter;                           //energy and pressures in the system.
double ACCEPTANCE_RATIO[2] = {0, 0};                 //0 - not accepted steps of rotation, 1 - accepted steps of rotation
int BALANCE_STEPS = 100;                             //steps for balance statistics
double delta = 0.5;                                  //MC parameter
double delta_angle = 90.0*(3.141592653589/180.0);     //MC parameter. Maximal rotation in rad
double R = 8.3144598;
double eps = 0.502e-21;                              // LJ energy for nitrogen in J
double N_a = 6.02214e+23;
double sigma = 331.8e-12;
double k_B = 1.38e-23;
double Rc = 5;                                             // Cut-off radius in sigma
double Rc2 = Rc*Rc;
double Qn2 = -4.453e-40/(sigma*sigma);             // Quadrupole moment of N2 molecule
const double eps0 = 8.85418781762e-12;                     // The permittivity of free space in C2 m-2 N-1
const double A = 1.0/(4.0*3.1415926535*eps0)/(sigma*sigma);    // Coulomb's constant
double C_q = A*(3.0/4.0)*pow(Qn2,2);
double dn2 = 0.33092224232;               // Distance between nitrogen atoms in sigma units
double dq1 = 0.25527426160;               // Distance between "+" charge and center of quadrupole in sigma units
double dq2 = 0.31464737794;               // Distance between "-" charge and center of quadrupole in sigma units
const double qe = 1.6021766208e-19;                   // The charge of one electron in C
double q = 0.373*qe;                            // Charge of the quadrupole points in C
double gm = 50;

//#include "writeConfigPBC.h"
#include "energies_and_forces.h"
//#include "initConfig.h"
#include "PBC2D.h"
#include "initConfigHerringbone.h"
//#include "initConfigPinwheel.h"
#include "PotentialEnergy.h"
//#include "Rosenbluth_algorithm_simple.h"
//#include "replace_the_trialParticle_and_update_energies.h"
#include "Metropolis_iteration.h"
#include "pressure_balance.h"
#include "layer_map.h"
#include "write_xy_matrix.h"
#include "writeData.h"
#include "write_xyz_file.h"

int main()
{
clock_t begin_time = clock();
 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 // Set configuration parameters
 double density = 0;                // Density
 double write_rad = 0.1098/0.3318;   // Disp. radius

 // N-N parameters: eps/k = 36.4 K, sigma = 0.3318 nm, Rc = 5 sigma
 // Temperature: 20 K
 // Reduced forms:
 // T* = kT/eps
 // F* = F*sigma/eps
 // P* = P*sigma3/eps

 // Set simulation parameters
 double temperature = 20;
 double Temp = temperature/36.4;               // Simulation temperature in units of eps/k
 double beta = 1.0/Temp;                       // Inverse temperature

 double Pt = 0;
 double press_X=0, press_Y=0, Energy=0;
 double en_2_av = 0;
 pressure press;
 press.X_LJ = 0;
 press.X_QQ = 0;
 press.Y_LJ = 0;
 press.Y_QQ = 0;

 // Set initial configuration
 double Lx=0,Ly=0;  // Linear size of the system
 int nPart = 400;
 vector <state> coordinates(nPart); // Vector of the molecules coordinates and angles

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Temperature" << "\t" << "Heat.Capacity" << "\t" << "E_per_molecule" << "\t" << "p_X" << "\t" << "p_Y" << "\t" << "p_X_LJ" << "\t" << "p_X_QQ" << "\t" << "p_Y_LJ" << "\t" << "p_Y_QQ" << "\t" << "Lx" << "\t" << "Ly" << endl;
 fileOutput.close();

 ////////////////////////////////////////////////////////////
 //         kMC simulation of systems with different N     //
 ////////////////////////////////////////////////////////////

 //for(int nPart = minPart; nPart < maxPart; nPart += stepPart)
 for(double state_dens = 10.5; state_dens < 10.6; state_dens += 1)
    {
     EN_AND_PR_counter.energy = 0;
     EN_AND_PR_counter.p.X_LJ = 0;
     EN_AND_PR_counter.p.X_QQ = 0;
     EN_AND_PR_counter.p.Y_LJ = 0;
     EN_AND_PR_counter.p.Y_QQ = 0;
     int frame = 0;
     // Generate an initial configuration for a fixed
     // number of particles and calculate required L

     //initConfig(nPart, density, sigma, coordinates, beta, Rc, A, C_q);   // Randomly distributed molecules
     initConfigHerringbone(nPart, density, coordinates, Lx, Ly, state_dens);       // Herringbone structure
     //initConfigPinwheel(nPart, density, coordinates, Lx, Ly, coeff);          // Pinwheel structure

     // Write the initial configuration
     write_xyz_file(nPart, Lx, Ly, temperature, coordinates, 0, 0, true); // Clear the xyz-file
     write_xyz_file(nPart, Lx, Ly, temperature, coordinates, frame, write_rad, false); // Write a current frame
     frame++;

     vector <vector <double>> xy_matrix(600, vector<double> (600));
     for(int i = 0; i < 600; i++){for(int j = 0; j < 600; j++){xy_matrix[i][j] = 0;}}

     // Calculate initial energy
     PotentialEnergy(nPart, Lx, Ly, coordinates, beta);


     /////////////////////////////
     // Set the Monte Carlo run //
     /////////////////////////////
     int nSteps = 5000;            // Total amount of MCS
     int nIter = nSteps * nPart;
     int nStepsEq = 2000;           // MCS for relaxation
     int nIterEq = nStepsEq * nPart;
     double Time = 0; // Total time of the equilibrium run
     double Mconf = 0; // Amount of configurations for chemical potential calculation with kMC
     double dt = 0;
     int balanceEq = 0;

     //////////////////////////////////////////////////
     //             Monte Carlo Simulation           //
     //////////////////////////////////////////////////

     for(int iter = 1; iter <= nIter; iter++)
        {
         if(((iter*100)%nIter)==0){cout << iter*100.0/nIter << " %" << endl;}

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
                press.X_LJ += EN_AND_PR_counter.p.X_LJ*dt;
                press.X_QQ += EN_AND_PR_counter.p.X_QQ*dt;
                press.Y_LJ += EN_AND_PR_counter.p.Y_LJ*dt;
                press.Y_QQ += EN_AND_PR_counter.p.Y_QQ*dt;
            if(((iter%(BALANCE_STEPS*nPart))==0 && iter != 0) || iter==nIterEq-1)
            {
                press.X_LJ /= Pt;
                press.X_QQ /= Pt;
                press.Y_LJ /= Pt;
                press.Y_QQ /= Pt;
                if (iter < 0.09*nIterEq && iter >= 0.05*nIterEq) { BALANCE_STEPS = 200; }
                if (iter < 0.15*nIterEq && iter >= 0.09*nIterEq) { BALANCE_STEPS = 300; }
                if (iter < 0.25*nIterEq && iter >= 0.15*nIterEq) { BALANCE_STEPS = 500; }
                if (iter < 0.46*nIterEq && iter >= 0.25*nIterEq) { BALANCE_STEPS = 1000; }
                if (iter >= 0.46*nIterEq) { BALANCE_STEPS = 2500; }

                pressure_balance ((press.X_LJ + press.X_QQ), (press.Y_LJ + press.Y_QQ), Lx, Ly, nPart, coordinates, beta);
                cout << float(clock ()-begin_time) /  CLOCKS_PER_SEC << endl;
                begin_time = clock();
                Pt = 0;
                press.X_LJ = 0;
                press.X_QQ = 0;
                press.Y_LJ = 0;
                press.Y_QQ = 0;
                balanceEq = 0;
                ACCEPTANCE_RATIO[0] = 0;
                ACCEPTANCE_RATIO[1] = 0;
            }
        }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

         // Collect the characteristics of interest at equilibrium
         if(iter > nIterEq)
           {


            if (iter == nIterEq+1)
            {
                Energy = 0;
                en_2_av = 0;
                press.X_LJ = 0;
                press.X_QQ = 0;
                press.Y_LJ = 0;
                press.Y_QQ = 0;
                Pt = 0;
            }

            //if (iter%(1000*nPart) == 0) {PotentialEnergy(nPart, Lx, Ly, coordinates, beta);}      //EN_AND_PR_counter correction
            // Update the total time and the amount of configurations for kMC (chemical potential calculation)
            if(rosenbluth) {Time += dt; Mconf++;}

            layer_map(nPart, coordinates, xy_matrix, Lx, Ly);
            //create xyz animation
            if((iter-1) % ((nIter-nIterEq)/300) == 0)
              {
                   write_xyz_file(nPart, Lx, Ly, temperature, coordinates, frame,  write_rad, false);
                   frame++;
              }


               Pt += dt;
               Energy += EN_AND_PR_counter.energy*dt;
               en_2_av += EN_AND_PR_counter.energy*EN_AND_PR_counter.energy*dt;
               press.X_LJ += EN_AND_PR_counter.p.X_LJ*dt;
               press.X_QQ += EN_AND_PR_counter.p.X_QQ*dt;
               press.Y_LJ += EN_AND_PR_counter.p.Y_LJ*dt;
               press.Y_QQ += EN_AND_PR_counter.p.Y_QQ*dt;
            }
         // A new random position is chosen uniformly over the whole volume of the system
         // and update the energies of all molecules in the system
         //if(rosenbluth) {replace_the_trialParticle_and_update_energies(nPart, trialPart, Lx, Ly, beta, coordinates);}
       }

     double mu = 0;
     if(rosenbluth) {mu = log(Mconf/Lx/Ly) - log(Time);}

            press.X_LJ = press.X_LJ/Pt;
            press.X_QQ = press.X_QQ/Pt;
            press.Y_LJ = press.Y_LJ/Pt;
            press.Y_QQ = press.Y_QQ/Pt;

            press.X_LJ = press.X_LJ/Lx/Ly/sigma/sigma*1000;
            press.X_QQ = press.X_QQ/Lx/Ly/sigma/sigma*1000;
            press.Y_LJ = press.Y_LJ/Lx/Ly/sigma/sigma*1000;
            press.Y_QQ = press.Y_QQ/Lx/Ly/sigma/sigma*1000;

            press_X = k_B*temperature*nPart/Ly/Lx/sigma/sigma*1000 + press.X_LJ + press.X_QQ;
            press_Y = k_B*temperature*nPart/Ly/Lx/sigma/sigma*1000 + press.Y_LJ + press.Y_QQ;

            Energy = Energy/Pt;
            en_2_av = en_2_av/Pt;

     // Write the final configuration
     //writeConfigPBC(nPart, density, sigma, Lx, Ly, coordinates, write_rad, "final");

     // Write the calculated data to a file
     writeData(temperature, (en_2_av-pow(Energy,2))/nPart, Energy/nPart*k_B*temperature*N_a/1000.0, press_X, press_Y, press.X_LJ, press.X_QQ, press.Y_LJ, press.Y_QQ, Lx, Ly);

     // Write the xy-matrix
     write_xy_matrix(nPart, Lx, Ly, temperature, xy_matrix);

     cout << "rho: " << density << "\t" << "mu: " << mu << "\t" << "en(kJ/mol): " << Energy/nPart*k_B*temperature*N_a/1000.0 << endl;

    }
 return 0;
}
