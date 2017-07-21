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
void set_state(double, double, double, double);
};

void state::set_state (double c_x, double c_y, double c_tetta, double c_phi) {
  x = c_x;
  y = c_y;
  tetta = c_tetta;
  phi = c_phi;
}

#include "distPBC.h"
#include "writeConfigPBC.h"
#include "Inter_potential.h"
#include "initConfig.h"
#include "PotentialEnergy.h"
#include "PBC2D.h"
#include "PotentialEnergyChange.h"

int main()
{

 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 // Set configuration parameters
 int nPart = 400;           // Number of N2 molecules
 double density = 0;        // Density
 double sigma = 1;          // LJ Diameter of nitrogen atom
 double write_rad = 0.1098;          // Disp. radius

 // N-N parameters: eps/k = 36.4 K, sigma = 0.3318 nm, Rc = 5 sigma
 // Temperature: 20 K
 // Reduced forms:
 // T* = kT/eps
 // F* = F*sigma/eps
 // P* = P*sigma3/eps

 // Set simulation parameters
 double Temp = 0.54945054945054945054945054945055;      // Simulation temperature in units of eps/k
 double beta = 1.0/Temp;                                // Inverse temperature
 double delta = 0.5;                                    // Maximal displacement in LJ units
 double delta_angle = 5;                                // Maximal rotation in degrees
 double Rc = 5;                                         // Cut-off radius in sigma
 double Qn2 = -4.453e-40;                               // Quadrupole moment of N2 molecule
 const double eps0 = 8.85418781762e-12;                 // The permittivity of free space in C2 m-2 N-1
 const double A = 1.0/(4.0*3.1415926535*eps0)/4.0;      // Coulomb's constant /4 to meet LJ calculation
 double C_q=A*3/4*pow(Qn2,2);                           // Coefficient of QQ interaction

 int nSteps = 10000;            // Total amount of MCS !!!
 int nIter = nSteps * nPart;
 int nStepsEq = 2000;          // MCS for relaxation
 int nIterEq = nStepsEq * nPart;

 // Write the model parameters to data-file
 stringstream name;
 name << "model_parameters" << ".dat";
 ofstream fileOutput(name.str().c_str(), ios_base::app);
 fileOutput << "Number of particles: " << nPart << endl;
 fileOutput << "Diameter of particles: " << sigma << endl;
 fileOutput << "T: " << Temp << endl;
 fileOutput.close();

 // Set initial configuration
 double Lx,Ly;  // Linear size of the system
 vector <state> coordinates(5000); // Vector of the molecules coordinates and angles

 for(density = 1.4; density < 1.41; density += 0.1)
    {

     // Generate an initial configuraton for a fixed
     // number of particles and calculate required L
     Lx = initConfig(nPart, density, sigma, coordinates, beta, Rc, A, C_q);
     Ly = Lx*2/sqrt(3);
     // Write the initial configuration
     writeConfigPBC(nPart, sigma, Lx, coordinates, write_rad, "initial");


     // Calculate initial energy
     double energy;
     energy = PotentialEnergy(nPart, sigma, Lx, Ly, Rc, coordinates, A, C_q);
     cout << "E: " << energy << endl;

     //////////////////////////////////////////////////
     //             Monte Carlo Simulation           //
     //////////////////////////////////////////////////

     //double vir = 0; // Lennard-Jones energy virial
     //double accept = 0; // Counter for accepatance rate
     int step = 0;
     //double w_sum = 0; // Widom's sum over exp(-dU/kT)
     //int mu_step = 0; // Counter for chemical potential calculation
     //int av_step = 0; // Counter for calculating averages

     for(int iter = 1; iter <= nIter; iter++)
        {
         step += 1;
         // Suggest a trial move for particle "trialPart"
         int trialPart;
         trialPart = RanGen.IRandom(0,(nPart-1));
         state trial_mol = coordinates[trialPart]; // Make a clone of trail particle

         if(step%2 == 0) // Move or Rotate a molecule
           {
            trial_mol.x = coordinates[trialPart].x + (2 * delta * RanGen.Random() - delta); // random(-delta; delta)
            trial_mol.y = coordinates[trialPart].y + (2 * delta * RanGen.Random() - delta);

            // Apply periodic boundary conditions
            trial_mol.x = PBC2D(Lx, trial_mol.x);
            trial_mol.y = PBC2D(Ly, trial_mol.y);

            // Calculate the change in energy due to this trial move
            double deltaE;
            deltaE = PotentialEnergyChange(nPart, sigma, Lx, Ly, Rc, coordinates, trial_mol, trialPart, A, C_q);

            if(RanGen.Random() < exp(-beta*deltaE))
              {
               // Accept displacement move
               coordinates[trialPart] = trial_mol;     // Update coordinates
               energy = energy + deltaE;   // Update energy
               //accept = accept + 1;
              }
           }
           else
           {
            trial_mol.tetta = coordinates[trialPart].tetta + (2 * delta_angle * RanGen.Random() - delta_angle); // random(-delta; delta)
            trial_mol.phi = coordinates[trialPart].phi + (2 * delta_angle * RanGen.Random() - delta_angle);

            // Calculate the change in energy due to this trial move
            double deltaE;
            deltaE = PotentialEnergyChange(nPart, sigma, Lx, Ly, Rc, coordinates, trial_mol, trialPart, A, C_q);

            if(RanGen.Random() < exp(-beta*deltaE))
              {
               // Accept displacement move
               coordinates[trialPart] = trial_mol;     // Update coordinates
               energy = energy + deltaE;   // Update energy
               //accept = accept + 1;
              }
           }
        }

     // Write the final configuration
     writeConfigPBC(nPart, sigma, Lx, coordinates, write_rad, "final");

         /*
         if(iter > nIterEq)
           {
            // Generate the coordinates of a test particle
            double test_x, test_y, test_z;
            test_x = L * RanGen.Random();
            test_y = L * RanGen.Random();
            test_z = L * RanGen.Random();
            w += LJ_Widom(test_x, test_y, test_z, L, nPart, beta, Rc, x, y, z);
            mu_step += 1;

            // Calculate the presuure virial
            if((iter%nPart) == 0)
              {
               av_step += 1;
               vir += LJ_EnergyVirial(nPart, L, Rc, x, y, z);
              }
           }

         // Tune the acceptance rate
         if(iter <= nIterEq && (iter%(100*nPart)) == 0){LJAcceptRatioControl(step, accept, delta);}
        }

     // Calculate the virial pressure
     vir = vir/av_step;             // Averaging the pressure
     double P = density*Temp + vir/3.0/pow(L, 3);

     // Claculate the chemical potential with Widom's method
     double mu = - log(w/mu_step)/beta;

     // Write the calculated data to a file
     writeData(nPart, L, "Density", "Pressure", "Chem. Potential", density, P, mu);

     // Write the final configuration
     writeConfigPBC(nPart, sigma, L, x, y, z, "final");

     // Calculate the acceptance rate
     accept = accept/step;

     cout << "rho: " << density << "\t" << "P: " << P << "\t" << "AR: " << accept << "\t" << "dr: " << delta << "\t" << "mu: " << mu << endl; */
    }

 return 0;
}
