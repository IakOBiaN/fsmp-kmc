// Random number generator.
// By default the seed is randomized from the wall clock (production runs).
// Define FSMP_RANDOM_SEED at compile time for deterministic, reproducible runs,
// e.g. for tests:  clang++ -DFSMP_RANDOM_SEED=12345 ...
#ifndef FSMP_RANDOM_SEED
#define FSMP_RANDOM_SEED ((int)time(0))
#endif
int seed = (FSMP_RANDOM_SEED);
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

class state {
public:
double x;
double y;
double phi;
double sin_phi;
double cos_phi;
double damping_coeff;
results ex_field_coeff;
double stat_weight;
results en_and_pr;
};

////////////// CONSTANTS /////////////////////////////////////////////////////

double R = 8.31446261815324;										// Gas constant in J per mol
double N_a = 6.02214076e+23;										//	Avogadro constant
double E_INF = 75.0;											// in kT units
const double PI = 3.14159265358979323846;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// GLOBAL VARIABLES ////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double temperature = temp_from;
double u_m = um_from;												// Parameter of the external field
results EN_AND_PR_counter;											// energy and pressures in the system.
double nPart_in_central_cell = 0;               					// molecules in central cell
double nPart_in_gas = 0;											// molecules in gas phase
double nPart_in_transition_zone = 0;								// molecules in nPart_in_transition_zone
double ACCEPTANCE_RATIO_r[2] = {0, 0};								// 0 - not accepted steps of rotation, 1 - accepted steps of rotation
double ACCEPTANCE_RATIO_m[2] = {0, 0};								// 0 - not accepted steps of move, 1 - accepted steps of move
int BALANCE_STEPS;													// steps for balance statistics
double density, gas_density, transition_zone_density;				// Actual density of the layer in mkMol per m^2
double lambda0 = sqrt(temperature / temperature_in_transition_zone);

bool HC_radius = false;                         					// Are we inside hard core radius (min_dist)?
bool findTrialPart = true;                      					// Condition for additional calculation of trialPart in kMC
int trialPart;
double sigma_2; // sigma in nm^2

// Forcefield: flat grid (distance, angle1, angle2) with ff_nang angular points per
// axis. Filled by read_forcefield; indexed through FF(i, j, k).
vector <double> forcefield;
int ff_nang;
double ff_fold_deg;
inline double & FF(int i, int j, int k)
{
	return forcefield[((size_t)i * ff_nang + j) * ff_nang + k];
}
// Minimal (hard core distance) and maximal distance between the molecules
double min_dist;
double min_dist_2;
double max_dist;
double max_dist_2;
int cut_index;
double dr;															// Delta between neighbor distances in the forcefield in A
double da;															// Delta between orientation angle of the single molecule
int frame = 0; // For visualization purpose

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "read_forcefield.h"
#include "unique_output_name.h"
#include "molecule_area.h"
#include "PBC2D.h"
#include "fields.h"
#include "stabilization_mask.h"
#include "interpolation.h"
#include "energies_and_forces_numerical.h"
#include "molecule_model.h"
#include "write_xyz_file.h"
#include "PotentialEnergy.h"
#include "Metropolis_iteration.h"
#include "Rosenbluth_iteration.h"
#include "pressure_balance.h"
#include "Weighted_averages.h"
#include "Widom_test.h"
#include "center_of_mass.h"
#include "StructureGenerator.h"

int main(int argc, char ** argv)
{
	if (argc == 2 && string(argv[1]) == "--version")
	{
		cout << "FSMP-kMC " << FSMP_VERSION << endl;
		return 0;
	}
	if (argc != 2)
	{
		cerr << "usage: " << argv[0] << " <parameters.txt>  (or --version)" << endl;
		cerr << "See configs/ for ready-to-run parameter files." << endl;
		return 2;
	}
	read_parameters(argv[1]);

	// Late initialization of the globals that depend on the parameters
	temperature = temp_from;
	u_m = um_from;
	lambda0 = sqrt(temperature / temperature_in_transition_zone);
	potential_name = p_name.c_str();
	if (param_seed_given) { RanGen.RandomInit(param_seed); }

	// The atomistic molecule model used for all xyz visualizations
	read_molecule_model(molecule_model_file);

	// The stabilization mask is damped by lambda(x) and must vanish in the gas
	// phase, otherwise the ideal-gas route to the chemical potential breaks
	if (stabilization_mask && lambda0 >= 1.0 && lambdam >= 1.0)
	{
		cerr << "ERROR: stabilization_mask = true requires a damping field (lambda < 1 somewhere): "
		     << "set temperature_in_transition_zone above the temperature or lambdam below 1" << endl;
		return 1;
	}

  complex_names();

	// Never overwrite the outputs of a previous run: shift the names if taken
	xyz_name = unique_output_name(xyz_name);
	if (structure_name == "calculate") { unit_cell_name = unique_output_name(unit_cell_name); }
	name_of_file_for_statistics.str(unique_output_name(name_of_file_for_statistics.str()));

 ///////////////////////////////////////
 //           Initialization          //
 ///////////////////////////////////////

 //////////////////////////////////////////////////////////////////////////////////////
 //////////////////////////// READ FORCEFIELD FROM FILE ///////////////////////////////
 /////////////////////////////////////////////////////////////////////////////////////

	// Read the precalculated numerical potential (binary; see tools/pack_forcefield).
	cout << "Now I'm reading the forcefield file." << endl;
	read_forcefield (potential_name, forcefield, ff_nang, min_dist, max_dist, dr, da, ff_fold_deg);

	min_dist_2 = min_dist*min_dist;
	max_dist_2 = max_dist*max_dist;
	cut_index = (int)(((max_dist - min_dist) / dr) + 0.5);

	// Reference area for the chemical potential (see sigma_mode in the parameter file)
	if (sigma_mode == "manual")
	{
		sigma_2 = sigma_manual * sigma_manual / 100.0;
	}
	else if (sigma_mode == "min_dist")
	{
		sigma_2 = min_dist_2 * PI / 4.0 / 100.0;
	}
	else // molecule_area
	{
		double area = molecule_area();
		sigma_2 = area / 100.0;
		cout << "Molecular area from the potential: " << area << " A^2" << endl;
	}
	cout << "Reference area sigma_2 (" << sigma_mode << "): " << sigma_2 << " nm^2" << endl;

 // Set configuration parameters

 double press_X = 0, press_Y = 0, Energy = 0, density = 0;
 double press_X_gas = 0, press_Y_gas = 0, Energy_gas = 0, gas_density = 0;
 double press_X_transition_zone = 0, Energy_transition_zone = 0, transition_zone_density = 0;
 double delta_p_over_interface = 0;
 double AR_r, AR_m;
 int N_test;																			// Counter for attempts to insert the test particle in Widom's algorythm
 double e_test;																	// Counter for energy change due to the insertion of the test particle

 ////////////////////////////////////////////////////////////
 //         MC simulation of systems with different N      //
 ////////////////////////////////////////////////////////////

 double Lx, Ly;  // Linear size of the system in A
 vector <state> coordinates(5000); // Vector of the molecules coordinates, angles and charges

	// Generating the initial structure for sequential MC simulation
  if (structure_name != "calculate")
  {
	   generate_structure(unit_cell_params, structure_name, coordinates, Lx, Ly);
  }
  else
  {
    calculate_unit_cell_params();
    generate_structure(unit_cell_params, coordinates, Lx, Ly);
    // The optimized cell has been reported and its animation written to
    // unit_cell_name; with optimize_only there is nothing left to do.
    if (optimize_only)
    {
      cout << endl << "optimize_only = true: stopping after the unit cell optimization." << endl;
      return 0;
    }
  }
  int nPart = unit_cell_params[0] * uc_in_x * uc_in_y;

	/////////////////////////////
	// Set the Monte Carlo run //
	/////////////////////////////
	int nIter = nSteps * nPart; // convert the MCS to iterations
	int nIterEq = nStepsEq * nPart;
	//vector <double> pressure_stat(nIter - nIterEq);
	//vector <double> energy_stat(nIter - nIterEq);

	// Clear up the xyz file
	write_xyz_file (xyz_name, nPart, density, Lx, Ly, temperature, coordinates, 0, 1, true);
	frame = 1;

	// Write the model parameters to data-file
	ofstream fileOutput(name_of_file_for_statistics.str().c_str(), ios_base::trunc);

	fileOutput << "Number of particles: " << nPart << endl;
	fileOutput << "Total number of MCS: " << nSteps << "  MCS for relaxation: " << nStepsEq << endl;
	fileOutput << "Maximal displacement, A: " << delta << "  Maximal rotation angle, deg: " << delta_angle << endl;
	fileOutput << "Lambda0: " << lambda0 << " Lambdam: " << lambdam << endl;

	fileOutput << "////////////////////////////////////////////////////////////////////////////////////////////" << endl << endl;

	// Columns are grouped by importance: the run conditions, the crystal state,
	// the chemical potential block, the gas phase with the analytical pressure,
	// and only then the strongly fluctuating virial pressures.
	fileOutput << "T, K" << "\t" << "u_m, kJ/mol" << "\t" << "lambda0" << "\t" << "lambda_m" << "\t"
	<< "Density, mkmol/m2" << "\t" << "Lx, A" << "\t" << "Ly, A" << "\t" << "Energy, kJ/mol" << "\t"
	<< "Excess chemical potential (ideal gas + u_m), kJ/mol" << "\t" << "Excess chemical potential (kMC, gas phase), kJ/mol" << "\t" << "Ideal excess part RTln(rho_gas*sigma2), kJ/mol" << "\t" << "Residual chemical potential (Widom), kJ/mol" << "\t" << "Excess chemical potential (ideal gas + Widom), kJ/mol" << "\t"
	<< "Gas density, mkmol/m2" << "\t" << "Analytical pressure in the crystal (Pg + dP), mN/m" << "\t" << "Pressure change over gas-solid interface (dP), mN/m" << "\t"
	<< "Virial total pressure, mN/m" << "\t" << "Virial excess pressure, mN/m" << "\t" << "Virial excess pressure along x, mN/m" << "\t" << "Virial excess pressure along y, mN/m" << endl;
	fileOutput.close();

bool u_m_loop_flag = true;
u_m = um_from;
um_step = abs(um_step);

 while (u_m_loop_flag)
 {

  bool temperature_loop_flag = true;
  temperature = temp_from;
  temp_step = abs(temp_step);

 while (temperature_loop_flag)
 {
	double beta = 1.0 / (R*temperature);  // Inverse temperature in units of (k_B*T)^-1
	// The transition-zone damping follows the current temperature; the stored
	// per-molecule coefficients are refreshed by PotentialEnergy below.
	lambda0 = sqrt(temperature/temperature_in_transition_zone);

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////// SYSTEM COUNTERS //////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	Energy = 0;
	press_X = 0;
	press_Y = 0;
	Energy_gas = 0;
	press_X_gas = 0;
	press_Y_gas = 0;
	Energy_transition_zone = 0;
	press_X_transition_zone = 0;
	delta_p_over_interface = 0;
	gas_density = 0;
	ACCEPTANCE_RATIO_r[0] = 0;  // rejected rotations
	ACCEPTANCE_RATIO_r[1] = 0;  // accepted rotations
	ACCEPTANCE_RATIO_m[0] = 0;  // rejected translations
	ACCEPTANCE_RATIO_m[1] = 0; // accepted translations
	int balanceEq = 0;
  double center_of_mass_x = 0;
////////////////////////////////////////////////////////////////////////////////////////////////////

// Calculate initial energy
	PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
  center_of_mass_x = center_of_mass(nPart, coordinates);
  weighted_averages_in_central_cell(coordinates, nPart, Lx, Ly);
	cout << endl << "_________INITIAL DATA_________" << endl;
	cout << endl << "u_m: " << u_m << endl;
	cout << endl << "Central cell" << endl;
	weighted_averages_in_central_cell(coordinates, nPart, Lx, Ly);
	cout << "Density: " << nPart_in_central_cell*(1.0e+26)/(Lx/4.0*Ly)/N_a << "\t" << " Energy: " << EN_AND_PR_counter.energy/1000.0/nPart_in_central_cell << "\t" << " P: " << (R*temperature*(1.0e+23)*nPart_in_central_cell/(Lx/4.0*Ly)/N_a)+((EN_AND_PR_counter.p_X + EN_AND_PR_counter.p_Y)/2.0/(Lx/4.0)/Ly*1e23/N_a)<< endl;
	cout << "P_X: " << (EN_AND_PR_counter.p_X/(Lx/4.0)/Ly*1e23/N_a) << "\t" << "P_Y: " << (EN_AND_PR_counter.p_Y/(Lx/4.0)/Ly*1e23/N_a) <<  endl;
	cout << endl;
	//////////////////////////////////////////////////
	//             Monte Carlo Simulation           //
	//////////////////////////////////////////////////
	double sum_iterations = 0;
	double percent = 0;
	double dt = 0; // Time step in kMC. For MMC it should always be 1.0
	double Pt = 0; // Integral time of the kinetic Monte Carlo run
	for(int iter = 1; iter <= nIter; iter++)
		{
			percent += 1;
			if(percent > nIter / 100.0)
				{
					frame++;
					write_xyz_file (xyz_name, nPart, density, Lx, Ly, temperature, coordinates, frame, 1, false);
					PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
					cout << int(iter*100.0/nIter) << " %" << endl;
					percent = 0;
				}

        //if it is begining (we need it for MMC too)
        if (iter < nIter * 0.02)
        {
          Metropolis_iteration(nPart, Lx, Ly, beta, coordinates, true, center_of_mass_x);
          dt = 1.0;
        }
        else
        {
          if (kMC)
          {
            //sometimes we need some MMC
            if (iter % 1000 == 0)
            {
              //some MMC with small steps
              for (int mmc_iter = 0; mmc_iter < 250; mmc_iter++)
              {
                Metropolis_iteration(nPart, Lx, Ly, beta, coordinates, true, center_of_mass_x);
              }
              //we don't have dt for kMC now
              findTrialPart = true;
            }
            //usual kMC iteration
            dt = Rosenbluth_iteration(Lx, Ly, nPart, coordinates, dt, beta, iter, trialPart, findTrialPart, center_of_mass_x);
            //after kMC iteration we have dt
            findTrialPart = false;
          }
          else
          {
              //two types of MMC steps (but we will average everything)
              if (iter % 10 == 0)
              {
                //usual MMC step
                Metropolis_iteration(nPart, Lx, Ly, beta, coordinates, true, center_of_mass_x);
              }
              else
              {
                //MMC step like in kMC
                Metropolis_iteration(nPart, Lx, Ly, beta, coordinates, false, center_of_mass_x);
              }
          }
        }

				balanceEq++;
				BALANCE_STEPS = 1000;
				if((iter < nIterEq) && (balanceEq > nPart*0.1*BALANCE_STEPS))
					{
						Pt += dt;
						sum_iterations += 1;
						weighted_averages_in_central_cell(coordinates, nPart, Lx, Ly);
						Energy += EN_AND_PR_counter.energy*dt;
						press_X += EN_AND_PR_counter.p_X*dt;
						press_Y += EN_AND_PR_counter.p_Y*dt;

            pressure_change_over_interface(coordinates, nPart, Lx, Ly);
            delta_p_over_interface += EN_AND_PR_counter.p_X*dt;

						weighted_averages_in_gas(coordinates, nPart, Lx, Ly);
						gas_density += nPart_in_gas*dt;

						if(((iter%(BALANCE_STEPS*nPart))==0 && iter != 0) || iter==nIterEq-1)
							{
								Energy /= Pt;
								press_X /= Pt;
								press_Y /= Pt;
								delta_p_over_interface /= Pt;
								gas_density /= Pt;
/*
								if (iter < 0.09*nIterEq && iter >= 0.05*nIterEq) { BALANCE_STEPS = 200; }
								if (iter < 0.15*nIterEq && iter >= 0.09*nIterEq) { BALANCE_STEPS = 300; }
								if (iter < 0.25*nIterEq && iter >= 0.15*nIterEq) { BALANCE_STEPS = 500; }
								if (iter < 0.46*nIterEq && iter >= 0.25*nIterEq) { BALANCE_STEPS = 1000; }
								if (iter >= 0.46*nIterEq) { BALANCE_STEPS = 2500; }
*/
								pressure_balance_ratio(Energy, press_X, press_Y, Lx, Ly, nPart, coordinates, beta);
//								pressure_balance_ratio_analytical(Energy, press_X, delta_p_over_interface, gas_density, Lx, Ly, nPart, coordinates, beta);
                if (constant_pressure)
                {
                  um_tunning_to_constant_pressure(Energy, u_m, delta_p_over_interface, Lx, Ly, nPart, coordinates, beta);
                }
								AR_r = ACCEPTANCE_RATIO_r[1]/(ACCEPTANCE_RATIO_r[0]+ACCEPTANCE_RATIO_r[1]);
								AR_m = ACCEPTANCE_RATIO_m[1]/(ACCEPTANCE_RATIO_m[0]+ACCEPTANCE_RATIO_m[1]);
                cout << "P_X_vir: " << press_X*(1.0/(Lx/4.0)/Ly*1e23)/N_a << "\t" << "P_Y_vir: " << press_Y*(1.0/(Lx/4.0)/Ly*1e23)/N_a <<  endl;
                cout << "P_X_an: " << (delta_p_over_interface*1e23/Ly/N_a) + R*temperature*gas_density/1000.0 << " Gas phase pressure: " << R*temperature*gas_density/1000.0 << endl;
								cout << "AR_m: " << AR_m << " delta: " << delta << " AR_r: " << AR_r << " delta_ang: " << delta_angle << endl;
								if (AR_r < 0.25 && delta_angle > 5.0)
								{delta_angle -= 1.0;}
								if (AR_r > 0.3 && delta_angle < 120.0)
								{delta_angle += 1.0;}
								if (AR_m < 0.25 && delta > 0.05)
								{delta -= 0.02;}
								if (AR_m > 0.3 && delta < 0.8)
								{delta += 0.02;}

								Pt = 0;
								sum_iterations = 0;
								Energy = 0;
								press_X = 0;
								press_Y = 0;
								delta_p_over_interface = 0;
								gas_density = 0;
								balanceEq = 0;
								ACCEPTANCE_RATIO_r[0] = 0;
								ACCEPTANCE_RATIO_r[1] = 0;
								ACCEPTANCE_RATIO_m[0] = 0;
								ACCEPTANCE_RATIO_m[1] = 0;
							}
					}

			if(iter > nIterEq)
				{
					if (iter == nIterEq+1)
						{
							density = 0; Energy = 0; press_X = 0; press_Y = 0; Pt = 0; sum_iterations = 0;
							gas_density = 0; Energy_gas = 0; press_X_gas = 0; press_Y_gas = 0;
							transition_zone_density = 0; Energy_transition_zone = 0; press_X_transition_zone = 0;
							delta_p_over_interface = 0;
							N_test = 0; e_test = 0;
						}

					Pt += dt;
          sum_iterations += 1;

					weighted_averages_in_central_cell(coordinates, nPart, Lx, Ly);
					density += nPart_in_central_cell*dt;
					Energy += EN_AND_PR_counter.energy*dt;
					press_X += EN_AND_PR_counter.p_X*dt;
					press_Y += EN_AND_PR_counter.p_Y*dt;
					//energy_stat[sum_iterations] = EN_AND_PR_counter.energy*dt;
					//pressure_stat[sum_iterations] = (EN_AND_PR_counter.p_X + EN_AND_PR_counter.p_Y)*dt/2.0;

					weighted_averages_in_gas(coordinates, nPart, Lx, Ly);
					gas_density += nPart_in_gas*dt;
					Energy_gas += EN_AND_PR_counter.energy*dt;
					press_X_gas += EN_AND_PR_counter.p_X*dt;
					press_Y_gas += EN_AND_PR_counter.p_Y*dt;

					weighted_averages_in_transition_zone(coordinates, nPart, Lx, Ly);
					transition_zone_density += nPart_in_transition_zone*dt;
					Energy_transition_zone += EN_AND_PR_counter.energy*dt;
					press_X_transition_zone += EN_AND_PR_counter.p_X*dt;

					pressure_change_over_interface(coordinates, nPart, Lx, Ly);
					delta_p_over_interface += EN_AND_PR_counter.p_X*dt;

					if (widom_test_index){Widom_test(nPart, coordinates, Lx, Ly, beta, N_test, e_test);}
				}
		}

	density /= Pt;
	Energy /= Pt;
	press_X /= Pt;
	press_Y /= Pt;
	density *= (1.0e+26)/((Lx/4.0)*Ly)/N_a;
	press_X *= (1.0/(Lx/4.0)/Ly*1e23)/N_a;
	press_Y *= (1.0/(Lx/4.0)/Ly*1e23)/N_a;
	gas_density /= Pt;
	Energy_gas /= Pt;
	press_X_gas /= Pt;
	press_Y_gas /= Pt;
	gas_density *= (1.0e+26)/((Lx/4.0)*Ly)/N_a;
	press_X_gas *= (1.0/(Lx/4.0)/Ly*1e23)/N_a;
	press_Y_gas *= (1.0/(Lx/4.0)/Ly*1e23)/N_a;
	transition_zone_density /= Pt;
	Energy_transition_zone /= Pt;
	press_X_transition_zone /= Pt;
	delta_p_over_interface /= Pt;
	transition_zone_density *= (1.0e+26)/((3.0*Lx/8.0)*Ly)/N_a;
	press_X_transition_zone *= (1.0/(3.0*Lx/16.0)/Ly*1e23)/N_a;
	delta_p_over_interface *= 1e23/Ly/N_a;
//	double mu_res_widom = log(N_test/(e_test))/beta/1000.0; // Residual chemical potential calculated by WTPI
	double mu_res_widom = 0;
	if (widom_test_index){mu_res_widom = log(N_test/(e_test))/beta/1000.0;} // Residual chemical potential calculated by WTPI
	double mu_ex_kMC = (log(sum_iterations/Lx/Ly) - log(Pt) + log(sigma_2 * 100))/beta/1000.0;

	cout << endl << "u_m: " << u_m << endl;
	cout << "Crystal Data" << endl;
	cout << "Density: " << density << " mkmol/m2 " << "Lx of central cell: " << Lx/4.0 << " Ly: " << Ly << endl;
	cout << "T: " << temperature << "K" << " Energy per molecule: " << Energy/1000.0/(density*Lx*Ly*N_a/4.0/1.0e+26) << " kJ/mol" << endl;
	cout << "Total pressure: " << R*temperature*density/1000.0 + (press_X + press_Y)/2.0 << " mN/m" << endl;
	cout << "P_ex_MC: " << (press_X + press_Y)/2.0 << " mN/m" << " P_ex_MC_X: " << press_X << " mN/m" << " P_ex_MC_Y: " << press_Y << " mN/m" << endl;
//	cout << "Chemical potential in the simulation cell: " << mu_ex + log(density*Lx*Ly*N_a/4.0/1.0e+26)/beta/1000.0 << " kJ/mol" << endl;

	cout << endl;
	cout << "Transition Zone Data" << endl;
	cout << "Transition zone density: " << transition_zone_density << " mikromol/m2" << " Energy per molecule in transition zone: " << Energy_transition_zone/1000.0/(transition_zone_density*3.0*Lx*Ly*N_a/8.0/1.0e+26) << " kJ/mol" << endl;
	cout << "Pressure along X in transition zone: " << R*temperature*transition_zone_density/1000.0 + press_X_transition_zone << " mN/m" << " Ideal gas impact in transition zone pressure: " << R*temperature*transition_zone_density/1000.0 << " mN/m" << endl;

	cout << endl;
	cout << "Gas Phase Data" << endl;
	cout << "Gas density: " << gas_density << " mikromol/m2" << " Gas energy per molecule: " << Energy_gas/1000.0/(gas_density*Lx*Ly*N_a/4.0/1.0e+26) << " kJ/mol" << endl;
	cout << "Gas pressure along X: " << R*temperature*gas_density/1000.0 + press_X_gas << " mN/m" << " Gas pressure along Y: " << R*temperature*gas_density/1000.0 + press_Y_gas << " mN/m" << endl;
	cout << "Widom's residual chemical potential in gas phase: " << mu_res_widom << " kJ/mol / " << "Widom's excess chemical potential in the gas phase: " << mu_res_widom + log(gas_density*N_a*1.0e-24*sigma_2)/beta/1000.0 << " kJ/mol" << endl;
	cout << "kMC's excess chemical potential in gas phase: " << mu_ex_kMC << " kJ/mol / " << endl;

	cout << endl;
	cout << "Pressure change over gas-solid interface (dP): " << delta_p_over_interface << " mN/m" << " Analytical pressure in the crystal (Pg + dP): " << R*temperature*gas_density/1000.0 + delta_p_over_interface << endl;

	double mu_id_gas = log(gas_density*N_a*1.0e-24*sigma_2)/beta/1000.0;   // ideal excess part, kJ/mol
	ofstream fileOutput(name_of_file_for_statistics.str().c_str(), ios_base::app);
	fileOutput  << temperature << "\t" << u_m/1000.0 << "\t" << lambda0 << "\t" << lambdam
	<< "\t" << density << "\t" << Lx << "\t" << Ly << "\t" << Energy/1000.0/(density*Lx*Ly*N_a/4.0/1.0e+26)
	<< "\t" << mu_id_gas + u_m/1000.0 << "\t" << mu_ex_kMC << "\t" << mu_id_gas << "\t" << mu_res_widom << "\t" << mu_id_gas + mu_res_widom
	<< "\t" << gas_density << "\t" << R*temperature*gas_density/1000.0 + delta_p_over_interface << "\t" << delta_p_over_interface
	<< "\t" << R*temperature*density/1000.0 + (press_X + press_Y)/2.0 << "\t" << (press_X + press_Y)/2.0 << "\t" << press_X << "\t" << press_Y << endl;
	fileOutput.close();


   if (temp_from < temp_to)
   {
     temperature += temp_step;
     if (temperature > temp_to)
     {
       temperature_loop_flag = false;
     }
   }
   else
   {
     temperature -= temp_step;
     if (temperature < temp_to)
     {
       temperature_loop_flag = false;
     }
   }
 // end of while loop for temperature
 }

   if (um_from < um_to)
   {
     u_m += um_step;
     if (u_m > um_to)
     {
       u_m_loop_flag = false;
     }
   }
   else
   {
     u_m -= um_step;
     if (u_m < um_to)
     {
       u_m_loop_flag = false;
     }
   }
// end of while loop for u_m
}

 return 0;
}
