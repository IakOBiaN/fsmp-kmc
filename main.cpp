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
#include <valarray>
using namespace std;

// Random number generator
int seed = (int)time(0);
CRandomSFMT0 RanGen(seed);

// Class contains the function descrabing the state of the molecule:
// 1) coordinates
// 2) angles
class state {
public:
double x;
double y;
double phi;
double energy;
double mob;
void set_state(double, double, double, double, double);
};

void state::set_state (double c_x, double c_y, double c_phi, double c_en, double c_mob) {
  x = c_x;
  y = c_y;
  phi = c_phi;
  energy = c_en;
  mob = c_mob;
}

double R = 8.3144598;
double eps = 0.515e-21;                         // LJ energy for nitrogen in J
double N_a = 6.02214e+23;
double sigma = 331.8e-12;
double k_B = 1.38e-23;
double Rc = 5;                                             // Cut-off radius in sigma
double Rc2 = Rc*Rc;
double Qn2 = -4.453e-40/(sigma*sigma);             // Quadrupole moment of N2 molecule
const double eps0 = 8.85418781762e-12;                     // The permittivity of free space in C2 m-2 N-1
const double A = 1.0/(4.0*3.1415926535*eps0)/(sigma*sigma);    // Coulomb's constant
double C_q = A*(3.0/4.0)*pow(Qn2,2);

//#include "writeConfigPBC.h"
#include "Inter_potential.h"
//#include "initConfig.h"
#include "PBC2D.h"
#include "initConfigHerringbone.h"
#include "initConfigPinwheel.h"
#include "PotentialEnergy.h"
#include "PotentialEnergy_SEP.h"
#include "Rosenbluth_algorithm_simple.h"
#include "replace_the_trialParticle_and_update_energies.h"
#include "Metropolis_iteration.h"
//#include "inter_for_pressure.h"
//#include "inter_for_pressure_ab.h"
#include "virial_pressure.h"
#include "pressure_balance.h"
#include "layer_map.h"
#include "write_xy_matrix.h"
#include "writeData.h"
#include "write_xyz_file.h"

int main()
{

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
 double press_X_LJ = 0, press_X_QQ = 0, press_Y_LJ = 0, press_Y_QQ = 0, press_X=0, press_Y=0, Energy=0;
 double p_X_LJ=0, p_X_QQ=0, p_Y_LJ=0, p_Y_QQ=0, E_per_Part=0;

 // Set initial configuration
 double Lx=0,Ly=0;  // Linear size of the system
 vector <state> coordinates(5000); // Vector of the molecules coordinates and angles

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Density" << "\t" << "Chem. Potential" << "\t" << "E per molecule" << "\t" << "p_X" << "\t" << "p_Y" << "\t" << "p_X_LJ" << "\t" << "p_X_QQ" << "\t" << "p_Y_LJ" << "\t" << "p_Y_QQ" << "\t" << "Lx" << "\t" << "Ly" << endl;
 fileOutput.close();

 ////////////////////////////////////////////////////////////
 //         kMC simulation of systems with different N     //
 ////////////////////////////////////////////////////////////

 //for(int nPart = minPart; nPart < maxPart; nPart += stepPart)
 int nPart = 400;
 for(double state_dens = 10.5; state_dens < 10.6; state_dens += 1)
    {
     bool rosenbluth = false;    // If rosenbluth = false then Metropolis algorithm works

     int frame = 0;
     // Generate an initial configuration for a fixed
     // number of particles and calculate required L

     //initConfig(nPart, density, sigma, coordinates, beta, Rc, A, C_q);   // Randomly distributed molecules
     initConfigHerringbone(nPart, density, coordinates, Lx, Ly, state_dens);       // Herringbone structure
     //initConfigPinwheel(nPart, density, coordinates, Lx, Ly, coeff);          // Pinwheel structure

     // Write the initial configuration
     //writeConfigPBC(nPart, density, sigma, Lx, Ly, coordinates, write_rad, "initial");
     write_xyz_file(nPart, Lx, Ly, temperature, coordinates, 0, 0, true);
     write_xyz_file(nPart, Lx, Ly, temperature, coordinates, frame,  write_rad, false);
     frame++;

     vector <vector <double>> xy_matrix(600, vector<double> (600));
     for(int i = 0; i < 600; i++){for(int j = 0; j < 600; j++){xy_matrix[i][j] = 0;}}

     // Calculate initial energy
     PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
     //for(int i = 0; i < nPart; i++){cout << "[" << i << "]: " << coordinates[i].energy << endl;}

     // Set the Monte Carlo run
     int nSteps = 20000;            // Total amount of MCS
     int nIter = nSteps * nPart;
     int nStepsEq = 10000;          // MCS for relaxation
     int nIterEq = nStepsEq * nPart;
     double Time = 0; // Total time of the equilibrium run
     double Mconf = 0; // Amount of configurations for chemical potential calculation with kMC
     double dt = 0;

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

         if(rosenbluth) {trialPart = Rosenbluth_algorithm_simple(nPart, coordinates, dt);}  // kMC trial particle and dt calculation
         else {Metropolis_iteration(nPart, Lx, Ly, beta, coordinates);dt = 1.0;}              // Make a MC iteration


         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Pressure balance

        if(iter < nIterEq && (iter%nPart) == 0)
        {
                p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0;
                virial_pressure(nPart, Lx, Ly, beta, coordinates, p_X_LJ, p_X_QQ, p_Y_LJ, p_Y_QQ);
                Pt += dt;
                press_X_LJ += p_X_LJ*dt;
                press_X_QQ += p_X_QQ*dt;
                press_Y_LJ += p_Y_LJ*dt;
                press_Y_QQ += p_Y_QQ*dt;

            if((iter%(500*nPart))==0 && iter != 0)
            {
                press_X_LJ /= Pt;
                press_X_QQ /= Pt;
                press_Y_LJ /= Pt;
                press_Y_QQ /= Pt;

                //pressure_balance((press_X_LJ + press_X_QQ), (press_Y_LJ + press_Y_QQ), Lx, Ly, nPart, coordinates, beta);
                Pt = 0;
                press_X_LJ = 0;
                press_X_QQ = 0;
                press_Y_LJ = 0;
                press_Y_QQ = 0;
            }
        }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

         // Collect the characteristics of interest at equilibrium
         if(iter > nIterEq)
           {
            if (iter == nIterEq+1)
            {
                press_X_LJ = 0;
                press_X_QQ = 0;
                press_Y_LJ = 0;
                press_Y_QQ = 0;
                Pt = 0;
            }
            // Update the total time and the amount of configurations for kMC (chemical potential calculation)
            if(rosenbluth) {Time += dt; Mconf++;}

            layer_map(nPart, coordinates, xy_matrix, Lx, Ly);
            //create xyz animation
            if((iter-1) % ((nIter-nIterEq)/300) == 0)
              {
                   write_xyz_file(nPart, Lx, Ly, temperature, coordinates, frame,  write_rad, false);
                   frame++;
              }

            if((iter%nPart) == 0)
              {
               Pt += dt;

               // Calculate the pressure tensor
               p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0; E_per_Part = 0;
               virial_pressure(nPart, Lx, Ly, beta, coordinates, p_X_LJ, p_X_QQ, p_Y_LJ, p_Y_QQ);
               for (int i = 0;i < nPart; i++) {E_per_Part += coordinates[i].energy;}
               E_per_Part /= nPart;
               press_X_LJ += p_X_LJ*dt;
               press_X_QQ += p_X_QQ*dt;
               press_Y_LJ += p_Y_LJ*dt;
               press_Y_QQ += p_Y_QQ*dt;
               Energy += E_per_Part*dt;
              }
            }

         // A new random position is chosen uniformly over the whole volume of the system
         // and update the energies of all molecules in the system
         if(rosenbluth) {replace_the_trialParticle_and_update_energies(nPart, trialPart, Lx, Ly, beta, coordinates);}
       }

     double mu = 0;
     if(rosenbluth) {mu = log(Mconf/Lx/Ly) - log(Time);}

            press_X_LJ = press_X_LJ/Pt;
            press_X_QQ = press_X_QQ/Pt;
            press_Y_LJ = press_Y_LJ/Pt;
            press_Y_QQ = press_Y_QQ/Pt;

            press_X_LJ = press_X_LJ/Ly/sigma/sigma*1000;
            press_X_QQ = press_X_QQ/Ly/sigma/sigma*1000;
            press_Y_LJ = press_Y_LJ/Lx/sigma/sigma*1000;
            press_Y_QQ = press_Y_QQ/Lx/sigma/sigma*1000;

            press_X = k_B*temperature*nPart/Ly/Lx/sigma/sigma*1000 + press_X_LJ + press_X_QQ;
            press_Y = k_B*temperature*nPart/Ly/Lx/sigma/sigma*1000 + press_Y_LJ + press_Y_QQ;

            Energy = Energy/Pt;

     // Write the final configuration
     //writeConfigPBC(nPart, density, sigma, Lx, Ly, coordinates, write_rad, "final");

     // Write the calculated data to a file
     writeData(density, mu, Energy*k_B*temperature*N_a/1000.0, press_X, press_Y, press_X_LJ, press_X_QQ, press_Y_LJ, press_Y_QQ, Lx, Ly);

     // Write the xy-matrix
     write_xy_matrix(nPart, Lx, Ly, temperature, xy_matrix);

     cout << "rho: " << density << "\t" << "mu: " << mu << "\t" << "en(kJ/mol): " << Energy*k_B*temperature*N_a/1000.0 << endl;

    }
 return 0;
}
