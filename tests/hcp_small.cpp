// Regression test: the full engine on the small committed grid (tests/data/).
// Single (T, u_m) point, two MC steps; the checked value is the deterministic
// INITIAL energy of the HCP crystal, printed before any Monte Carlo move.
// The grid is a strided (coarse) copy of TMA_simple_2020, so the energy is a
// regression pin, not a physically converged number.
#include "../includes.h"

double temp_from = 300;
double temp_to = 300;
double temp_step = 10;
double um_from = 0.0;
double um_to = 0.0;
double um_step = 5000.0;
double temperature_in_transition_zone = 900;
double lambdam = 0.0;
int nSteps = 2;       // tiny: we only need the initial energy print
int nStepsEq = 1;
bool constant_pressure = false;
double constant_pressure_value = 0;
bool kMC = true;

string p_name = "data/TMA_simple_2020_s4.v2.bin";
string structure_name = "TMA_HCP_simple_2020";
int uc_in_x = 22;
int uc_in_y = 6;
double free_space = 0.24;
int total_molecule_directions = 3;
double angle_1 = 120;
double angle_2 = 240;
double delta = 2.0;
double delta_angle = 60.0;
bool widom_test_index = false;

string unit_cell_name = "0_hcp_small.xyz";
string xyz_name = "1_hcp_small.xyz";
stringstream name_of_file_for_statistics;
void complex_names()
{
  name_of_file_for_statistics << "2_hcp_small.dat";
}

vector<double> unit_cell_params;
void calculate_unit_cell_params()
{
  unit_cell_params.push_back(2);
  unit_cell_params.push_back(11.100);
  unit_cell_params.push_back(19.200);
  unit_cell_params.push_back(0); unit_cell_params.push_back(0); unit_cell_params.push_back(90);
  unit_cell_params.push_back(11.089); unit_cell_params.push_back(59.967); unit_cell_params.push_back(90);
}

const char * potential_name = p_name.c_str();

#include "../program_body.cpp"
