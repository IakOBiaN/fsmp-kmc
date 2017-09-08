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
double tetta;
double phi;
double energy;
double mob;
void set_state(double, double, double, double, double, double);
};

void state::set_state (double c_x, double c_y, double c_tetta, double c_phi, double c_en, double c_mob) {
  x = c_x;
  y = c_y;
  tetta = c_tetta;
  phi = c_phi;
  energy = c_en;
  mob = c_mob;
}

#include "writeConfigPBC.h"
#include "Inter_potential.h"
//#include "initConfig.h"
#include "PBC2D.h"
#include "initConfigHerringbone.h"
#include "initConfigPinwheel.h"
#include "PotentialEnergy.h"
#include "Rosenbluth_algorithm_simple.h"
#include "replace_the_trialParticle_and_update_energies.h"
#include "Metropolis_iteration.h"
//#include "inter_for_pressure.h"
#include "inter_for_pressure_ab.h"
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
 int minPart = 400;
 int maxPart = 420;
 int stepPart = 20;
 double density = 0;        // Density
                            // in units of the commensurate (C)
                            // in-plane herringbone structure,
                            // 0.0636 molecules N2 per A^2
 double write_rad = 0.1098/0.3318; // Disp. radius

 // N-N parameters: eps/k = 36.4 K, sigma = 0.3318 nm, Rc = 5 sigma
 // Temperature: 20 K
 // Reduced forms:
 // T* = kT/eps
 // F* = F*sigma/eps
 // P* = P*sigma3/eps

 // Set simulation parameters
 double temperature = 20;
 double Temp = temperature/36.4;                            // Simulation temperature in units of eps/k
 double beta = 1.0/Temp;                                    // Inverse temperature
 double Rc = 5;                                             // Cut-off radius in sigma
 double Rc2 = Rc*Rc;
 double Qn2 = -4.453e-40;                                   // Quadrupole moment of N2 molecule
 const double eps0 = 8.85418781762e-12;                     // The permittivity of free space in C2 m-2 N-1
 const double A = 1.0/(4.0*3.1415926535*eps0)/331.8e-12;    // Coulomb's constant
 double C_q=A*(3/4)*pow(Qn2,2);                             // Coefficient of QQ interaction

 double Pt = 0;
 double press_N = 0, press_T = 0, press = 0, Ener = 0;
 double p_N, p_T, p_Tot, E_per_Part;

 // Set initial configuration
 double Lx=0,Ly=0;  // Linear size of the system
 vector <state> coordinates(5000); // Vector of the molecules coordinates and angles

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Density" << "\t" << "Chem. Potential" << "\t" << "E per molecule" << "\t" << "p" << "\t" << "p_N" << "\t" << "p_T" << "\t" << "Lx" << "\t" << "Ly" << endl;
 fileOutput.close();

 ////////////////////////////////////////////////////////////
 //         kMC simulation of systems with different N     //
 ////////////////////////////////////////////////////////////

 //for(int nPart = minPart; nPart < maxPart; nPart += stepPart)
 int nPart = 144;
 for(double coeff = 1.06; coeff < 1.061; coeff += 0.01)
    {
     bool rosenbluth = true;    // If rosenbluth = false then Metropolis algorithm works

     int frame = 0;
     // Generate an initial configuration for a fixed
     // number of particles and calculate required L

     //initConfig(nPart, density, sigma, coordinates, beta, Rc, A, C_q);   // Randomly distributed molecules
     initConfigHerringbone(nPart, density, coordinates, Lx, Ly, coeff);       // Herringbone structure
     //initConfigPinwheel(nPart, density, coordinates, Lx, Ly, coeff);          // Pinwheel structure

     // Write the initial configuration
     //writeConfigPBC(nPart, density, sigma, Lx, Ly, coordinates, write_rad, "initial");
     write_xyz_file(nPart, Lx, Ly, temperature, coordinates, 0, 0, true);

     vector <vector <double>> xy_matrix(600, vector<double> (600));
     for(int i = 0; i < 600; i++){for(int j = 0; j < 600; j++){xy_matrix[i][j] = 0;}}

     // Calculate initial energy
     PotentialEnergy(nPart, Lx, Ly, Rc, Rc2, coordinates, A, C_q, beta);
     //for(int i = 0; i < nPart; i++){cout << "[" << i << "]: " << coordinates[i].energy << endl;}

     // Set the Monte Carlo run
     int nSteps = 50000;            // Total amount of MCS
     int nIter = nSteps * nPart;
     int nStepsEq = 10000;          // MCS for relaxation
     int nIterEq = nStepsEq * nPart;
     double Time = 0; // Total time of the equilibrium run
     double Mconf = 0; // Amount of configurations for chemical potential calculation with kMC
     double nMetroConf = 0; // Amount of equilibrium configurations in Metropolis run
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

         if(rosenbluth == true){trialPart = Rosenbluth_algorithm_simple(nPart, coordinates, dt);}
         // Make a Metropolis iteration
         else {Metropolis_iteration(nPart, Rc, Rc2, Lx, Ly, beta, A, C_q, coordinates);}


         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Pressure balance
/*
        if(iter < nIterEq && (iter%nPart) == 0)
        {
            Pt += dt;
            p_N = 0; p_T = 0;
            virial_pressure(nPart, Lx, Ly, beta, Rc, Rc2, coordinates, p_N, p_T, p_Tot, A, C_q);
            press_N += p_N*dt;
            press_T += p_T*dt;

            if((iter%(100*nPart))==0 && iter != 0)
            {
                press_N = press_N/Pt;                                 // Average normal pressure
                press_N = - press_N;
                press_T = press_T/Pt;                                 // Average normal pressure
                press_T = - press_T;
                pressure_balance(press_N, press_T, Lx, Ly, nPart, coordinates, Rc, Rc2, A, C_q, beta);
                Pt = 0;
                press_N = 0;
                press_T = 0;
            }
        }
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if((iter-1) % (nIter/300) == 0)
              {
                   write_xyz_file(nPart, Lx, Ly, temperature, coordinates, frame,  write_rad, false);
                   frame++;
              }

         // Collect the characteristics of interest at equilibrium
         if(iter > nIterEq)
           {
            // Update the total time
            if(rosenbluth == true){Time += dt;}
            // Update the amount of configurations
            Mconf += 1;

            layer_map(nPart, coordinates, xy_matrix);

            if((iter%nPart) == 0)
              {
               if(rosenbluth == true){Pt += dt;}

               // Calculate the pressure tensor
               p_N = 0; p_T = 0; p_Tot = 0; E_per_Part = 0;
               virial_pressure(nPart, Lx, Ly, beta, Rc, Rc2, coordinates, p_N, p_T, p_Tot, A, C_q);

               // Calculate the energy per molecule
               for(int i = 0; i < nPart; i++){E_per_Part += coordinates[i].energy;}
               E_per_Part /= nPart;

               if(rosenbluth == true)
               { // For calculating time averages with kMC
                press_N += p_N*dt;
                press_T += p_T*dt;
                press += p_Tot*dt;
                Ener += E_per_Part*dt;
               }
               else { // For calculating ensemble averages with sMC
                     press_N += p_N;
                     press_T += p_T;
                     press += p_Tot;
                     Ener += E_per_Part;
                     nMetroConf += 1;
                    }
              }
            }

         // A new random position is chosen uniformly over the whole volume of the system
         // and update the energies of all molecules in the system
         if(rosenbluth == true){replace_the_trialParticle_and_update_energies(nPart, trialPart, Rc, Rc2, Lx, Ly, beta, A, C_q, coordinates);}
       }

     double mu = 0;
     if(rosenbluth == true)
       {
         mu = log(Mconf/Lx/Ly) - log(Time);

         // kMC simulation
         // Calculate the virial pressure
         press_N = press_N/Pt;                   // Time average normal pressure
         press_N = Temp*(1.0 - press_N/density)/1000;
         press_T = press_T/Pt;                   // Time average tangential pressure
         press_T = Temp*(1.0 - press_T/density)/1000;
         press = press/Pt;                       // Time average total pressure
         press = Temp*(1.0 - press/2.0/density)/1000;
         Ener = Ener/Pt;
       }
       else { // Metropolis run
             press_N = press_N/nMetroConf;       // Ensemble average normal pressure
             press_N = Temp*(1.0 - press_N/density)/1000;
             press_T = press_T/nMetroConf;       // Ensemble average tangential pressure
             press_T = Temp*(1.0 - press_T/density)/1000;
             press = press/nMetroConf;           // Ensemble average total pressure
             press = Temp*(1.0 - press/2.0/density)/1000;
             Ener = Ener/nMetroConf;
            }

     // Write the final configuration
     //writeConfigPBC(nPart, density, sigma, Lx, Ly, coordinates, write_rad, "final");

     // Write the calculated data to a file
     writeData(density, mu, Ener, press, press_N, press_T, Lx, Ly);

     // Write the xy-matrix
     write_xy_matrix(nPart, Lx, Ly, temperature, xy_matrix);

     cout << "rho: " << density << "\t" << "mu: " << mu << "\t" << "e: " << Ener << endl;

    }
 return 0;
}
