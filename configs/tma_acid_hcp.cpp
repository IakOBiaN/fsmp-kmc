#include "../includes.h"
///////////////////////////////////////////////////////////////////////////////
//////////////// CONFIGURATION ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//MAIN PARAMETERS
//temperature in crystal (K)
double temp_from = 300;
double temp_to = 301;
double temp_step = 10;
//parameter of the external field
double um_from = 0.0;
double um_to = -100000.0;
double um_step = 5000.0;
//parameter of the external field
double temperature_in_transition_zone = 900;						// K
double lambdam = 0.0;
/////////// Set the lenght of MC run ///////////////////////////////
int nSteps = 500000;            // Total amount of MCS
int nStepsEq = 250000;           // MCS for relaxation
bool constant_pressure = false;
double constant_pressure_value = 0;
bool kMC = true;

//MODEL SETTINGS
//numerical potential name
string p_name = "../forcefields/TMA_simple_2020.v2.bin";
//you can use your own structures if set "structure_name" to "calculate"
string structure_name = "TMA_HCP_simple_2020";  //from StructureGenerator.h file
int uc_in_x = 22;  //unit cells in x direction
int uc_in_y = 6;   //unit cells in y direction
double free_space = 0.24;  //persent of the free space in the elongated cell
/////// molecules visualization block /////////////
int total_molecule_directions = 3;
double angle_1 = 120;
double angle_2 = 240;
//ADDITIONAL PARAMETERS
double delta = 2.0;													// MC parameter. Maximal shift of the molecule
double delta_angle = 60.0;    										// MC parameter. Maximal rotation in degrees
bool widom_test_index = false;

//name of xyz file of unit cell optimization (if you want to optimize a unit cell)
string unit_cell_name = "0_calculate_animation.xyz";
//name of xyz file for visualization
string xyz_name = "1_xyz_for_calculations_TMA_hcp.xyz";
stringstream name_of_file_for_statistics;
void complex_names()
{
  name_of_file_for_statistics << "2_statistics_" << structure_name << "_" << "Xuc" << uc_in_x << "_" << "Yuc" << uc_in_y << "_" << "FreeSp" << free_space << ".dat";
}

//IF YOU WANT TO CALCULATE UNIT CELL
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
