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
double um_step = 1000.0;
//parameter of the external field
double temperature_in_transition_zone = 900;						// K
double lambdam = 0.0;
/////////// Set the lenght of MC run ///////////////////////////////
int nSteps = 1000000;            // Total amount of MCS
int nStepsEq = 500000;           // MCS for relaxation
bool constant_pressure = false;
double constant_pressure_value = 0;
bool kMC = true;

//MODEL SETTINGS
//numerical potential name
string p_name = "../forcefields/TPA_qPBE_crystal_Dhb5.0_r5.5_16_dr01_da5.dat";
//you can use your own structures if set "structure_name" to "calculate"
string structure_name = "TPA_horizontal_ladder_qPBE_Dreiding_Dhb5.0";  //from StructureGenerator.h file
int uc_in_x = 30;  //unit cells in x direction
int uc_in_y = 7;   //unit cells in y direction
double free_space = 0.24;  //persent of the free space in the elongated cell
/////// molecules visualization block (now it is 2 directions) /////////////
int total_molecule_directions = 2;
//for TPA
double angle_1 = 180;
double angle_2 = 240;       // it can be any angle as long as total_molecule_directions is 2

//ADDITIONAL PARAMETERS
double delta = 2.0;													// MC parameter. Maximal shift of the molecule
double delta_angle = 60.0;    										// MC parameter. Maximal rotation in degrees
bool widom_test_index = false;

//name of xyz file of unit cell optimization (if you want to optimize a unit cell)
string unit_cell_name = "0_calculate_animation.xyz";
//name of xyz file for visualization
string xyz_name = "1_xyz_for_calculations_TPA_ladder.xyz";
stringstream name_of_file_for_statistics;
void complex_names()
{
  name_of_file_for_statistics << "2_statistics_" << structure_name << "_" << "Xuc" << uc_in_x << "_" << "Yuc" << uc_in_y << "_" << "FreeSp" << free_space << ".dat";
}

//IF YOU WANT TO CALCULATE UNIT CELL
vector<double> unit_cell_params;

void calculate_unit_cell_params()
{
  //molecules count
  unit_cell_params.push_back(2);
  //size along x axis for unit cell
  unit_cell_params.push_back(11.735);
  //size along y axis for unit cell
   unit_cell_params.push_back(9.647);

  //unit cell graph (graph distance, graph angle, molecule angle)
  //first molecule
  unit_cell_params.push_back(-0.000888541);
  unit_cell_params.push_back(90.137772);
  unit_cell_params.push_back(90.0);
  //second molecule
  unit_cell_params.push_back(7.303);
  unit_cell_params.push_back(143.002);
  unit_cell_params.push_back(90.253273);
}

const char * potential_name = p_name.c_str();

#include "../program_body.cpp"
