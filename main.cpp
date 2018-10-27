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
double energy_LJ;
double energy_QQ;
double p_X_LJ, p_Y_LJ, p_X_QQ, p_Y_QQ;
results();      //constructor
results operator+(const results& b) {
         results res;
         res.energy = this->energy + b.energy;
         res.energy_LJ = this->energy_LJ + b.energy_LJ;
         res.energy_QQ = this->energy_QQ + b.energy_QQ;
         res.p_X_LJ = this->p_X_LJ + b.p_X_LJ;
         res.p_Y_LJ = this->p_Y_LJ + b.p_Y_LJ;
         res.p_X_QQ = this->p_X_QQ + b.p_X_QQ;
         res.p_Y_QQ = this->p_Y_QQ + b.p_Y_QQ;
         return res;
      }
results operator-(const results& b) {
         results res;
         res.energy = this->energy - b.energy;
         res.energy_LJ = this->energy_LJ - b.energy_LJ;
         res.energy_QQ = this->energy_QQ - b.energy_QQ;
         res.p_X_LJ = this->p_X_LJ - b.p_X_LJ;
         res.p_Y_LJ = this->p_Y_LJ - b.p_Y_LJ;
         res.p_X_QQ = this->p_X_QQ - b.p_X_QQ;
         res.p_Y_QQ = this->p_Y_QQ - b.p_Y_QQ;
         return res;
      }
};

//constructor
results::results(void) {
   energy = 0;
   energy_LJ = 0;
   energy_QQ = 0;
   p_X_LJ = 0;
   p_Y_LJ = 0;
   p_X_QQ = 0;
   p_Y_QQ = 0;
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
results EN_AND_PR_counter;                           //energy and pressures in the system.
double ACCEPTANCE_RATIO_r[2] = {0, 0};               //0 - not accepted steps of rotation, 1 - accepted steps of rotation
double ACCEPTANCE_RATIO_m[2] = {0, 0};               //0 - not accepted steps of move, 1 - accepted steps of move
int BALANCE_STEPS = 100;                             //steps for balance statistics
double sigma = 3.318;
double delta = 0.3;                            //MC parameter. Maximal shift of the molecule
double delta_angle = 30.0;    //MC parameter. Maximal rotation in degrees
double R = 8.3144598;
double N_a = 6.02214e+23;
double k_B = 1.38e-23;
const double PI  =3.141592653589793238463;
//double gm = 50;
double density = 0;
double temperature = 0;

#include "read_forcefield.h"
// Forcefield for TMA-TMA pair
vector <vector <vector <double> > > forcefield;
vector <vector <vector <double> > > energy_LJ;
vector <vector <vector <double> > > energy_QQ;
// Minimal and maximal distance between the molecules (hard core distance)
double min_dist,max_dist;
// Delta between neighbor distances in the forcefield in A
double dr;
// Delta between orientation angle of the single molecule
double da;
int frame = 0;

#include "interpolation.h"
#include "energies_and_forces.h"
#include "energies_and_forces_2.h"
#include "energies_and_pressures.h"
#include "PBC2D.h"
#include "initConfigHerringbone.h"
#include "PotentialEnergy.h"
#include "PotentialEnergy_2.h"
#include "Metropolis_iteration.h"
#include "pressure_balance.h"
#include "layer_map.h"
#include "write_xy_matrix.h"
#include "writeData.h"
#include "write_xyz_file.h"

int main()
{

  // Fill in the forcefield
  // First dimension is distance
  // Second dimension is angle of first molecule
  // Third dimension is angle of second molecule
  for (int i = 0; i < 150; i++) {
      vector< vector<double> > mat; // Create an empty matrix
      for (int j = 0; j < 361; j++) {
          vector<double> row; // Create an empty row
              for (int k =0; k <361; k++) {
                  row.push_back(0);
              }
          mat.push_back(row); // Add an element (column) to the row
      }
      forcefield.push_back(mat); // Add the row to the main vector
      energy_LJ.push_back(mat);
      energy_QQ.push_back(mat);
  }
  // Read the forcefield from "forcefield.dat"
  cout << "read forcefield.dat file" << endl;
  read_forcefield ("forcefield.dat", forcefield, min_dist,max_dist, dr, da);
  cout << "read force_LJ.dat file" << endl;
  read_forcefield ("energy_LJ.dat", energy_LJ, min_dist,max_dist, dr, da);
  cout << "read force_QQ.dat file" << endl;
  read_forcefield ("energy_QQ.dat", energy_QQ, min_dist,max_dist, dr, da);
 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 // Set configuration parameters
 double density = 0;                // Density

 double Pt = 0;
 double press_X = 0, press_Y = 0, press_X_LJ = 0, press_Y_LJ = 0,press_X_QQ = 0, press_Y_QQ = 0, Energy=0;
 double en_2_av = 0;
 double cap_n = 1.0;
 double fluent_capacity = 0;
 double persent = 0;

 /////////////////////////////
 // Set the Monte Carlo run //
 /////////////////////////////
 int nPart = 100;
 int nSteps = 100000;            // Total amount of MCS
 int nIter = nSteps * nPart;
 int nStepsEq = 50000;           // MCS for relaxation
 int nIterEq = nStepsEq * nPart;
 double Lx=0,Ly=0;  // Linear size of the system
 double state_dens = 10.5; // mkMol of N2 per m^2
 vector <state> coordinates(nPart); // Vector of the molecules coordinates and angles

 // Write the model parameters to data-file
 stringstream name;
 name <<  "statistics.dat";
 ofstream fileOutput(name.str().c_str(), ios_base::trunc);
 fileOutput << "Temperature" << "\t" << "Heat.Capacity(reccurent)" << "\t" << "Heat.Capacity" << "\t" << "E_per_molecule" << "\t" << "p_X" << "\t" << "p_Y" << "\t" << "p_X_LJ" << "\t"<< "p_X_QQ" << "\t"<< "p_Y_LJ" << "\t"<< "p_Y_QQ" << "\t" << "Lx" << "\t" << "Ly" << endl;
 fileOutput.close();

 ////////////////////////////////////////////////////////////
 //         kMC simulation of systems with different N     //
 ////////////////////////////////////////////////////////////

 //for(int nPart = minPart; nPart < maxPart; nPart += stepPart)

for(temperature = 20; temperature < 31; temperature += 2.0)
    {
      //Generete a random distribution of TMA molecules at fixed density
     initConfigHerringbone(nPart, density, coordinates, Lx, Ly, state_dens);
     write_xyz_file_N2 (nPart, Lx, Ly, temperature, coordinates, 0, 1, true);
     EN_AND_PR_counter.energy = 0;
     EN_AND_PR_counter.energy_LJ = 0;
     EN_AND_PR_counter.energy_QQ = 0;
     EN_AND_PR_counter.p_X_LJ = 0;
     EN_AND_PR_counter.p_Y_LJ = 0;
     EN_AND_PR_counter.p_X_QQ = 0;
     EN_AND_PR_counter.p_Y_QQ = 0;
	 Pt = 0;
   press_X = 0;
   press_Y = 0;
	 press_X_LJ = 0;
	 press_Y_LJ = 0;
   press_X_QQ = 0;
	 press_Y_QQ = 0;
	 ACCEPTANCE_RATIO_r[0] = 0;
	 ACCEPTANCE_RATIO_r[1] = 0;
	 ACCEPTANCE_RATIO_m[0] = 0;
	 ACCEPTANCE_RATIO_m[1] = 0;
	 double Time = 0; // Total time of the equilibrium run
	 double Mconf = 0; // Amount of configurations for chemical potential calculation with kMC
	 double dt = 0;
	 int balanceEq = 0;
     persent = 0;

	 double beta = 1.0 / (k_B*temperature);  // Inverse temperature in units of 1/J
   int mol1_number = 0,mol2_number = 1;
   coordinates[mol1_number].x = 0.0;
   coordinates[mol1_number].y = 0.0;
   coordinates[mol2_number].x = 2.8055;
   coordinates[mol2_number].y = 1.855;
   results test_exact = energies_and_forces_2(coordinates[mol1_number], coordinates[mol2_number], Lx, Ly,beta);
   results test_approx = energies_and_forces(coordinates[mol1_number], coordinates[mol2_number], Lx, Ly,beta);
   //coordinates[mol2_number].x = 4.5;
   cout << "COORDINATES=" << coordinates[mol1_number].x << " " << coordinates[mol1_number].y << " two=" << coordinates[mol2_number].x << " " << coordinates[mol2_number].y << endl;
   cout << "exact_press=" << test_exact.p_X_LJ + test_exact.p_X_QQ << " " << test_exact.p_Y_LJ + test_exact.p_Y_QQ << endl;
   cout << "approx_press=" << test_approx.p_X_LJ + test_approx.p_X_QQ << " " << test_approx.p_Y_LJ + test_approx.p_Y_QQ << endl;

   double diff = 0.01;
   double test_p_X,test_p_Y;
   state mol1 = coordinates[mol1_number];
   state mol2 = coordinates[mol2_number];
   mol1.x = coordinates[mol1_number].x + diff;
   mol2.x = coordinates[mol2_number].x - diff;
   double energy_one = energies_and_forces(mol1, mol2, Lx, Ly,beta).energy;
   mol1.x = coordinates[mol1_number].x - diff;
   mol2.x = coordinates[mol2_number].x + diff;
   double energy_two = energies_and_forces(mol1, mol2, Lx, Ly,beta).energy;
   test_p_X = (energy_one - energy_two)/(diff*2.0)*(coordinates[mol2_number].x-coordinates[mol1_number].x)/2.0;

   mol1 = coordinates[mol1_number];
   mol2 = coordinates[mol2_number];
   mol1.y = coordinates[mol1_number].y + diff;
   mol2.y = coordinates[mol2_number].y - diff;
   energy_one = energies_and_forces(mol1, mol2, Lx, Ly,beta).energy;
   mol1.y = coordinates[mol1_number].y - diff;
   mol2.y = coordinates[mol2_number].y + diff;
   energy_two = energies_and_forces(mol1, mol2, Lx, Ly,beta).energy;
   test_p_Y = (energy_one - energy_two)/(diff*2.0)*(coordinates[mol2_number].y-coordinates[mol1_number].y)/2.0;

   cout << "new_pressure=" << test_p_X << " " << test_p_Y << endl;
   break;
     vector <vector <double> > xy_matrix(1000, vector<double> (1000));
     for(int i = 0; i < 1000; i++){for(int j = 0; j < 1000; j++){xy_matrix[i][j] = 0;}}
     // Calculate initial energy
     PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
     cout << "energy=" << (EN_AND_PR_counter.energy/1000.0)*(N_a/nPart)/beta << endl;
     //PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
     //cout << "approx_energy=" << EN_AND_PR_counter.energy << endl;

     //////////////////////////////////////////////////
     //             Monte Carlo Simulation           //
     //////////////////////////////////////////////////

     for(int iter = 1; iter <= nIter; iter++)
        {

          if((iter%(nPart*50) == 0) || (iter == 1))
          {
           frame++;
           write_xyz_file_N2 (nPart, Lx, Ly, temperature, coordinates, frame, 1.094, false);
          }
         persent += 1;
         if(persent > nIter/100.0)
         {
           double old_energy = EN_AND_PR_counter.energy;
           PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
           if (abs(old_energy - EN_AND_PR_counter.energy) > abs(0.01*EN_AND_PR_counter.energy))
           {
              cout << "AHTUNG!!! Energy problem. Exact=" << EN_AND_PR_counter.energy << " diff_energy=" << old_energy << endl;
           }
           cout << int(iter*100.0/nIter) << " %" << endl;persent = 0;
         }

         int trialPart;
         // Choose a particle to be displaced
         // according to the ROSENBLUTH scheme
         // and calculate the duration of the current configuration
         //PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
         if(!rosenbluth) {Metropolis_iteration(nPart, Lx, Ly, beta, coordinates); dt = 1.0;}     // Make a MC iteration
         //else {trialPart = Rosenbluth_algorithm_simple(nPart, coordinates, dt);}                 // kMC trial particle and dt calculation


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Pressure balance

        balanceEq++;
        if((iter < nIterEq) && (balanceEq > nPart*0.1*BALANCE_STEPS))
        {
            Pt += dt;
            press_X += EN_AND_PR_counter.p_X_LJ*dt + EN_AND_PR_counter.p_X_QQ*dt;
            press_Y += EN_AND_PR_counter.p_Y_LJ*dt + EN_AND_PR_counter.p_Y_QQ*dt;
            if(((iter%(BALANCE_STEPS*nPart))==0 && iter != 0) || iter==nIterEq-1)
            {
                press_X /= Pt;
                press_Y /= Pt;
                if (iter < 0.09*nIterEq && iter >= 0.05*nIterEq) { BALANCE_STEPS = 200; }
                if (iter < 0.15*nIterEq && iter >= 0.09*nIterEq) { BALANCE_STEPS = 300; }
                if (iter < 0.25*nIterEq && iter >= 0.15*nIterEq) { BALANCE_STEPS = 500; }
                if (iter < 0.46*nIterEq && iter >= 0.25*nIterEq) { BALANCE_STEPS = 1000; }
                if (iter >= 0.46*nIterEq) { BALANCE_STEPS = 2500; }

                //pressure_balance (press_X, press_Y, Lx, Ly, nPart, coordinates, beta);
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
                press_X_LJ = 0;
                press_Y_LJ = 0;
                press_X_QQ = 0;
                press_Y_QQ = 0;
                Pt = 0;
            }

            // Update the total time and the amount of configurations for kMC (chemical potential calculation)
            if(rosenbluth) {Time += dt; Mconf++;}

            layer_map_N2(nPart, coordinates, xy_matrix, Lx, Ly);
            if (cap_n > 1) {fluent_capacity = (cap_n - 1.0)/cap_n*fluent_capacity + (cap_n-1.0)/(cap_n*cap_n)*(Energy/(Pt+1) - EN_AND_PR_counter.energy)*(Energy/(Pt+1) - EN_AND_PR_counter.energy);}

            cap_n++;
            Pt += dt;
            Energy += EN_AND_PR_counter.energy*dt;
            en_2_av += EN_AND_PR_counter.energy*EN_AND_PR_counter.energy*dt;
            press_X_LJ += EN_AND_PR_counter.p_X_LJ*dt;
            press_Y_LJ += EN_AND_PR_counter.p_Y_LJ*dt;
            press_X_QQ += EN_AND_PR_counter.p_X_QQ*dt;
            press_Y_QQ += EN_AND_PR_counter.p_Y_QQ*dt;
           }
         // A new random position is chosen uniformly over the whole volume of the system
         // and update the energies of all molecules in the system
         //if(rosenbluth) {replace_the_trialParticle_and_update_energies(nPart, trialPart, Lx, Ly, beta, coordinates);}
       }

     double mu = 0;
     if(rosenbluth) {mu = log(Mconf/Lx/Ly) - log(Time);}

            press_X_LJ /= Pt;
            press_Y_LJ /= Pt;
            press_X_QQ /= Pt;
            press_Y_QQ /= Pt;

            press_X_LJ *= (1.0/Lx/Ly*1e20*1000)/beta;  //it means p_x_lj = p_x_lj/Lx/Ly/sigma/sigma*1e20*1000 mN/m
            press_Y_LJ *= (1.0/Lx/Ly*1e20*1000)/beta;
            press_X_QQ *= (1.0/Lx/Ly*1e20*1000)/beta;
            press_Y_QQ *= (1.0/Lx/Ly*1e20*1000)/beta;

            /*press_X = R*temperature*nPart/Ly/Lx + press_X_LJ + press_X_QQ;
            press_Y = R*temperature*nPart/Ly/Lx + press_Y_LJ + press_Y_QQ;*/

            Energy = Energy/Pt;
            en_2_av = en_2_av/Pt;

     // Write the calculated data to a file
     writeData(temperature, fluent_capacity/nPart, (en_2_av-pow(Energy,2))/nPart, (Energy/1000.0)*(N_a/nPart)/beta, press_X_LJ+press_X_QQ, press_Y_LJ+press_Y_QQ,press_X_LJ, press_X_QQ, press_Y_LJ,press_Y_QQ, Lx, Ly,0.0);

     // Write the xy-matrix
     write_xy_matrix(nPart, Lx, Ly, temperature, xy_matrix);

     cout << "rho: " << density << " mkMol/m^2 \t" << "mu: " << mu << "\t" << "en: " << (Energy/1000.0)*(N_a/nPart)/beta << " kJ/mol" << endl;
     cout << "temp=" << temperature << " P_X_LJ=" << press_X_LJ << " P_Y_LJ=" << press_Y_LJ << " P_X_QQ=" << press_X_QQ << " P_Y_QQ=" << press_Y_QQ << endl;

    }
 return 0;
}
