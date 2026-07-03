// FSMP-kMC entry point. Build once, run with a parameter file:
//
//   make
//   ./fsmp.out configs/tma_acid_hcp.txt
//
// All simulation parameters are read at run time by read_parameters.h; the
// files in configs/ are ready-to-run examples and document every key. Paths
// inside a parameter file are relative to the directory the program is started
// from, and all output files are written there too.
#include "includes.h"

// Simulation parameters (filled from the parameter file; the zero/empty values
// here are placeholders, the parser enforces that required keys are present)
double temp_from = 0, temp_to = 0, temp_step = 0;      // crystal temperature loop, K
double um_from = 0, um_to = 0, um_step = 0;            // external field loop, J/mol
double temperature_in_transition_zone = 0;             // K
double lambdam = 0;
int nSteps = 0;                                        // total MCS
int nStepsEq = 0;                                      // MCS for relaxation
bool constant_pressure = false;
double constant_pressure_value = 0;
bool kMC = true;
string p_name;                                         // numerical potential (binary v2)
string structure_name;                                 // named structure or "calculate"
int uc_in_x = 0, uc_in_y = 0;                          // unit cells along x and y
double free_space = 0;                                 // fraction of free space in the elongated cell
int total_molecule_directions = 0;                     // visualization
double angle_1 = 0, angle_2 = 0;                       // visualization
double delta = 0;                                      // maximal MC shift, A
double delta_angle = 0;                                // maximal MC rotation, deg
bool widom_test_index = false;
string unit_cell_name;                                 // xyz animation of the cell optimization
string xyz_name;                                       // xyz trajectory
string statistics_name;                                // optional; default is built below
vector<double> unit_cell_params;                       // filled by the unit_cell key (calculate)

#include "read_parameters.h"

stringstream name_of_file_for_statistics;
void complex_names()
{
	if (statistics_name.empty())
	{
		name_of_file_for_statistics << "2_statistics_" << structure_name << "_" << "Xuc" << uc_in_x
		                            << "_" << "Yuc" << uc_in_y << "_" << "FreeSp" << free_space << ".dat";
	}
	else
	{
		name_of_file_for_statistics << statistics_name;
	}
}

void calculate_unit_cell_params()
{
	// unit_cell_params is already filled by the parser (the unit_cell key);
	// its presence and shape were validated there.
}

const char * potential_name = 0;   // set in main once p_name is known

#include "program_body.cpp"
