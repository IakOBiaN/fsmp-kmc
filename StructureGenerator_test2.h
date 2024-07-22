using namespace std;

void generate_elongated_cell(vector <double> &params, vector <state> &coordinates, double &Lx, double &Ly)
{
  results empty_field;
  double x_uc = params[1];
  double y_uc = params[2];

  double full_Lx = x_uc * uc_in_x;
  Lx = full_Lx / (1 - 2.0 * free_space);
  double shift_of_structure = (Lx - full_Lx) / 2.0;
  Ly = y_uc * uc_in_y;
  int molecules = 0;

  for(int i = 0; i < uc_in_x; i++)
  {
    for(int j = 0; j < uc_in_y; j++)
    {
      coordinates[molecules].x = shift_of_structure + i * x_uc + params[3] * cos(params[4] / 180.0 * PI);
      coordinates[molecules].y = j * y_uc + params[3] * sin(params[4] / 180.0 * PI);
      coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
      coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
      coordinates[molecules].phi = params[5];
      coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
      coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
      molecules++;
      for (int mol = 1; mol < params[0]; mol++)
      {
        int number = 3 + mol * 3;
        coordinates[molecules].x = coordinates[molecules - 1].x + params[number] * cos(params[number + 1] / 180.0 * PI);
        coordinates[molecules].y = coordinates[molecules - 1].y + params[number] * sin(params[number + 1] / 180.0 * PI);
        coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
        coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
        coordinates[molecules].phi = params[number + 2];
        coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
        coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
        molecules++;
      }
    }
  }

  double mass_center = center_of_mass(molecules, coordinates);
  shift_of_structure = Lx / 2.0 - mass_center;

  for (int i = 0; i < molecules; i++)
  {
    coordinates[i].x += shift_of_structure;
    coordinates[i].damping_coeff = damping_field(coordinates[i].x, Lx); // Lambda^1/2
    coordinates[i].ex_field_coeff = external_field(coordinates[i].x, Lx); // u_ext
    coordinates[i].stat_weight = weights_for_central_cell (coordinates[i].x, Lx);
  }

  cout << endl << "Elongated cell was generated: " << endl;
  cout << "N: " << molecules << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}

void generate_structure(vector <double> &params, string structure_name, vector <state> &coordinates, double &Lx, double &Ly)
{
  //old optimization
	if (structure_name == "TMA_fCW_qB3LYP_PBE_Dreiding_Dhb5.4")
	{
    unit_cell_params.push_back(6);
    unit_cell_params.push_back(17.2785);
    unit_cell_params.push_back(29.858);
    unit_cell_params.push_back(-1.43758);
    unit_cell_params.push_back(0.294197);
    unit_cell_params.push_back(33.1213);
    unit_cell_params.push_back(9.95879);
    unit_cell_params.push_back(29.927);
    unit_cell_params.push_back(93.1209);
    unit_cell_params.push_back(9.92876);
    unit_cell_params.push_back(149.515);
    unit_cell_params.push_back(97.0724);
    unit_cell_params.push_back(9.9157);
    unit_cell_params.push_back(29.5732);
    unit_cell_params.push_back(33.6068);
    unit_cell_params.push_back(9.91344);
    unit_cell_params.push_back(149.785);
    unit_cell_params.push_back(93.7853);
    unit_cell_params.push_back(9.42189);
    unit_cell_params.push_back(20.2558);
    unit_cell_params.push_back(90.3681);
	}
  //optimized
  if (structure_name == "TMA_CW_qB3LYP_PBE_Dreiding_Dhb5.4")
  {
     unit_cell_params.push_back(4);
     unit_cell_params.push_back(17.258);
     unit_cell_params.push_back(29.948);

     unit_cell_params.push_back(0.235903);
     unit_cell_params.push_back(-0.764381);
     unit_cell_params.push_back(27.0098);

     unit_cell_params.push_back(9.96516);
     unit_cell_params.push_back(30.0011);
     unit_cell_params.push_back(87.0006);

     unit_cell_params.push_back(9.96484);
     unit_cell_params.push_back(90.1917);
     unit_cell_params.push_back(27.192);

     unit_cell_params.push_back(9.99148);
     unit_cell_params.push_back(149.835);
     unit_cell_params.push_back(87.5855);
  }
  //optimized
  if (structure_name == "TMA_SF_qB3LYP_PBE_Dreiding_Dhb5.4")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.76494);
     unit_cell_params.push_back(16.9659);
     unit_cell_params.push_back(-1.4425);
     unit_cell_params.push_back(1.15426);
     unit_cell_params.push_back(38.7484);
     unit_cell_params.push_back(9.77999);
     unit_cell_params.push_back(60.1552);
     unit_cell_params.push_back(39.2479);
  }
  //optimized
  if (structure_name == "TPA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.90997);
     unit_cell_params.push_back(11.698);
     unit_cell_params.push_back(0.771123);
     unit_cell_params.push_back(-0.837901);
     unit_cell_params.push_back(1.00002);
     unit_cell_params.push_back(7.59232);
     unit_cell_params.push_back(50.475);
     unit_cell_params.push_back(1.00003);
  }
  //optimized
  if (structure_name == "TPA_vertical_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(11.6946);
     unit_cell_params.push_back(9.90999);
     unit_cell_params.push_back(-0.74105);
     unit_cell_params.push_back(91.1906);
     unit_cell_params.push_back(91);
     unit_cell_params.push_back(7.58993);
     unit_cell_params.push_back(140.459);
     unit_cell_params.push_back(90.9999);
  }
  if (structure_name == "TPA_horizontal_ladder_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.85863);
     unit_cell_params.push_back(14.9832);
     unit_cell_params.push_back(0.460148);
     unit_cell_params.push_back(2.47113);
     unit_cell_params.push_back(2.56578);
     unit_cell_params.push_back(9.41467);
     unit_cell_params.push_back(51.4863);
     unit_cell_params.push_back(80.2502);
  }
  if (structure_name == "TPA_vertical_ladder_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(15.0487);
     unit_cell_params.push_back(9.85162);
     unit_cell_params.push_back(3.16403);
     unit_cell_params.push_back(0.749759);
     unit_cell_params.push_back(93.0274);
     unit_cell_params.push_back(9.65954);
     unit_cell_params.push_back(37.2725);
     unit_cell_params.push_back(-11.2231);
  }
  //optimized
  if (structure_name == "IPA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(16.9632);
     unit_cell_params.push_back(7.33756);
     unit_cell_params.push_back(3.1684);
     unit_cell_params.push_back(-2.30584);
     unit_cell_params.push_back(147.618);
     unit_cell_params.push_back(9.85106);
     unit_cell_params.push_back(30.6955);
     unit_cell_params.push_back(327.696);
  }
	//optimized
	if (structure_name == "IPA_vertical_chain_qPBE_Dreiding_Dhb5.0")
	{
		 unit_cell_params.push_back(2);
		 unit_cell_params.push_back(7.41172);
		 unit_cell_params.push_back(16.7164);
		 unit_cell_params.push_back(0.472696);
		 unit_cell_params.push_back(-0.646095);
		 unit_cell_params.push_back(59.7553);
		 unit_cell_params.push_back(9.84645);
		 unit_cell_params.push_back(57.8756);
		 unit_cell_params.push_back(240.191);
	}
	// optimized
	if (structure_name == "IPA_hexagones_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(24);
		unit_cell_params.push_back(42.6379);
		unit_cell_params.push_back(49.1998);
		//1
		unit_cell_params.push_back(11.4791);
		unit_cell_params.push_back(0.0714106);
		unit_cell_params.push_back(238.214);
		//2
		unit_cell_params.push_back(9.88775);
		unit_cell_params.push_back(120.151);
		unit_cell_params.push_back(298.021);
		//3
		unit_cell_params.push_back(7.69802);
		unit_cell_params.push_back(89.9078);
		unit_cell_params.push_back(178.57);
		//4
		unit_cell_params.push_back(9.87203);
		unit_cell_params.push_back(59.6657);
		unit_cell_params.push_back(237.711);
		//5
		unit_cell_params.push_back(9.87488);
		unit_cell_params.push_back(120.171);
		unit_cell_params.push_back(299.081);
		//6
		unit_cell_params.push_back(7.33438);
		unit_cell_params.push_back(89.9816);
		unit_cell_params.push_back(177.82);
		//7
		unit_cell_params.push_back(12.41);
		unit_cell_params.push_back(22.2247);
		unit_cell_params.push_back(-2.32758);
		//8
		unit_cell_params.push_back(9.86637);
		unit_cell_params.push_back(240.134);
		unit_cell_params.push_back(58.0192);
		//9
		unit_cell_params.push_back(9.90949);
		unit_cell_params.push_back(300.02);
		unit_cell_params.push_back(118.865);
		//10
		unit_cell_params.push_back(7.55);
		unit_cell_params.push_back(270.095);
		unit_cell_params.push_back(-1.97112);
		//11
		unit_cell_params.push_back(9.85007);
		unit_cell_params.push_back(240.122);
		unit_cell_params.push_back(57.8898);
		//12
		unit_cell_params.push_back(9.85007);
		unit_cell_params.push_back(299.999);
		unit_cell_params.push_back(118.172);
		//13
		unit_cell_params.push_back(9.85);
		unit_cell_params.push_back(0.659447);
		unit_cell_params.push_back(177.63);
		//14
		unit_cell_params.push_back(9.86994);
		unit_cell_params.push_back(59.9726);
		unit_cell_params.push_back(238.041);
		//15
		unit_cell_params.push_back(9.8696);
		unit_cell_params.push_back(120);
		unit_cell_params.push_back(298.249);
		//16
		unit_cell_params.push_back(7.55138);
		unit_cell_params.push_back(90.0516);
		unit_cell_params.push_back(178.538);
		//17
		unit_cell_params.push_back(9.85537);
		unit_cell_params.push_back(60.0484);
		unit_cell_params.push_back(238.541);
		//18
		unit_cell_params.push_back(9.93114);
		unit_cell_params.push_back(120.218);
		unit_cell_params.push_back(298.191);
		//19
		unit_cell_params.push_back(7.54902);
		unit_cell_params.push_back(29.177);
		unit_cell_params.push_back(57.7393);
		//20
		unit_cell_params.push_back(9.85);
		unit_cell_params.push_back(299.651);
		unit_cell_params.push_back(117.65);
		//21
		unit_cell_params.push_back(7.52084);
		unit_cell_params.push_back(269.726);
		unit_cell_params.push_back(-1.65575);
		//22
		unit_cell_params.push_back(9.85);
		unit_cell_params.push_back(239.844);
		unit_cell_params.push_back(57.5787);
		//23
		unit_cell_params.push_back(9.85);
		unit_cell_params.push_back(300.191);
		unit_cell_params.push_back(117.997);
		//24
		unit_cell_params.push_back(7.50352);
		unit_cell_params.push_back(269.753);
		unit_cell_params.push_back(-1.95989);
	}
	//optimized
	if (structure_name == "PA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(2);
	  unit_cell_params.push_back(12.48);
	  unit_cell_params.push_back(9.85);
		unit_cell_params.push_back(0.547938);
		unit_cell_params.push_back(0.0439552);
		unit_cell_params.push_back(127.03);
		unit_cell_params.push_back(9.93052);
		unit_cell_params.push_back(51.4783);
		unit_cell_params.push_back(307.031);
	}
	//optimized
	if (structure_name == "PA_vertical_chain_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(2);
		unit_cell_params.push_back(9.43172);
		unit_cell_params.push_back(12.8256);
		unit_cell_params.push_back(2.95425);
		unit_cell_params.push_back(0.146893);
		unit_cell_params.push_back(35.2501);
		unit_cell_params.push_back(9.92587);
		unit_cell_params.push_back(40.2124);
		unit_cell_params.push_back(214.654);
	}

  generate_elongated_cell(params, coordinates, Lx, Ly);
}

void generate_structure(vector <double> &params, vector <state> &coordinates, double &Lx, double &Ly)
{
  double temp_E_INF = E_INF;
  E_INF = 1e200;
  results empty_field;
  results en_and_press;
  int params_amount = params.size() - 1;
  bool first = true;
	bool unit_cell_size_variation = false;
  double beta = 1.0 / (R * 300);
  int N = params[0] * 3 * 3;
  double id_P = 0, Press = 0, Px = 0, Py = 0;
	double temp_energy = 0;

	//filling the unit cell for the initial unit cell parameters
	Lx = params[1] * 3;
	Ly = params[2] * 3;
	double x_uc = params[1];
	double y_uc = params[2];

	/////// OPTIMIZATION PARAMENERS //////////////////////
	double delta_dist_and_angle = (x_uc + y_uc)*0.1/2.0;
	double delta_uc = (x_uc + y_uc)*0.01/2.0;
	//////////////////////////////////////////////////////

	int molecules = 0;
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			coordinates[molecules].x = i * x_uc + params[3] * cos(params[4] / 180.0 * PI);
			coordinates[molecules].y = j * y_uc + params[3] * sin(params[4] / 180.0 * PI);
			coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
			coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
			coordinates[molecules].phi = params[5];
			coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
			coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
			coordinates[molecules].damping_coeff = 1.0;
			coordinates[molecules].ex_field_coeff = empty_field;
			coordinates[molecules].stat_weight = 1.0;
			coordinates[molecules].en_and_pr = empty_field;
			molecules++;
			for (int mol = 1; mol < params[0]; mol++)
			{
				int number = 3 + mol * 3;
				coordinates[molecules].x = coordinates[molecules - 1].x + params[number] * cos(params[number + 1] / 180.0 * PI);
				coordinates[molecules].y = coordinates[molecules - 1].y + params[number] * sin(params[number + 1] / 180.0 * PI);
				coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
				coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
				coordinates[molecules].phi = params[number + 2];
				coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
				coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
				coordinates[molecules].damping_coeff = 1.0;
				coordinates[molecules].ex_field_coeff = empty_field;
				coordinates[molecules].stat_weight = 1.0;
				coordinates[molecules].en_and_pr = empty_field;
				molecules++;
			}
		}
	}

	write_xyz_file (unit_cell_name, N, density, Lx, Ly, temperature, coordinates, 0, 1, 1);

	// Calculate the energy of the initial unit cell
	for(int molA = 0; molA < (N - 1); molA++)
	{
		for(int molB = (molA + 1); molB < N; molB++)
			{
				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
				en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
				coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
				coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
			}
	}

	// claculate properties initial unit cell
	weighted_averages_in_central_cell(coordinates, N, Lx, Ly);

	// Assert on initial structure
	if (HC_radius){cout << "ERROR!!! HC_radius!!!" << endl;}

	Px = (EN_AND_PR_counter.p_X/(Lx/4.0)/Ly*1e23/N_a);
	Py = (EN_AND_PR_counter.p_Y/(Lx/4.0)/Ly*1e23/N_a);
	id_P = R*temperature*nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a/1000.0;
	Press = id_P + (Px + Py)/2.0;
	temp_energy = EN_AND_PR_counter.energy;

	cout << endl << "\\\\\\\\" << "Initial properties" << "\\\\\\\\" << endl;
	cout << " Density" << "\t" << " Energy" << "\t" << "P" << "\t" << " P_X" << "\t" << "P_Y" <<  endl;
	cout << nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a << "\t" << " Energy: " << EN_AND_PR_counter.energy / 1000.0 / nPart_in_central_cell << "\t" << Press << "\t" << Px << "\t" << Py <<  endl << endl;
	cout << "Step" << "\t" << "Iterations" << "\t" << "UC_change" << "\t" << "Shift_or_rotation" << "\t"
					<< "Lx" << "\t" << "Ly" << "\t"
					<< " Density" << "\t" << " Energy" << "\t" << "P" << "\t" << " P_X" << "\t" << "P_Y" <<  endl;

	int step = 0;
  int counter = 0;

  while (counter < 100000)
	{
		HC_radius = false;
		unit_cell_size_variation = false;
		bool energy_decision = false;
    int param_number;
    double trial_delta_dist_and_angle = 0;
		double trial_delta_uc = 0;\

		// Read the current parameters
		Lx = params[1] * 3;
		Ly = params[2] * 3;
		x_uc = params[1];
		y_uc = params[2];
		double default_Lx = Lx;
		double default_Ly = Ly;
		double default_x_uc = x_uc;
		double default_y_uc = y_uc;


		int molecules = 0;
    //filling the unit cell for the current unit cell parameters
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        coordinates[molecules].x = i * x_uc + params[3] * cos(params[4] / 180.0 * PI);
        coordinates[molecules].y = j * y_uc + params[3] * sin(params[4] / 180.0 * PI);
        coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
        coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
        coordinates[molecules].phi = params[5];
        coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
        coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
        coordinates[molecules].damping_coeff = 1.0;
        coordinates[molecules].ex_field_coeff = empty_field;
        coordinates[molecules].stat_weight = 1.0;
        coordinates[molecules].en_and_pr = empty_field;
        molecules++;
        for (int mol = 1; mol < params[0]; mol++)
        {
          int number = 3 + mol * 3;
          coordinates[molecules].x = coordinates[molecules - 1].x + params[number] * cos(params[number + 1] / 180.0 * PI);
          coordinates[molecules].y = coordinates[molecules - 1].y + params[number] * sin(params[number + 1] / 180.0 * PI);
          coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
          coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
          coordinates[molecules].phi = params[number + 2];
          coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
          coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
          coordinates[molecules].damping_coeff = 1.0;
          coordinates[molecules].ex_field_coeff = empty_field;
          coordinates[molecules].stat_weight = 1.0;
          coordinates[molecules].en_and_pr = empty_field;
          molecules++;
        }
      }
    }
		// Write the current state if it was accepted on previous step of the algorithm
		if(counter == 0) {write_xyz_file (unit_cell_name, N, density, Lx, Ly, temperature, coordinates, 0, 1, 0);}

		// Calculate the energy of the current state
		for(int molA = 0; molA < (N - 1); molA++)
		{
			for(int molB = (molA + 1); molB < N; molB++)
				{
					en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
					en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
					coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
					coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
				}
		}
		weighted_averages_in_central_cell(coordinates, N, Lx, Ly);
		temp_energy = EN_AND_PR_counter.energy; // energy of the current state

		// Calculate the pressures for the current state
		Px = (EN_AND_PR_counter.p_X/(Lx/4.0)/Ly*1e23/N_a);
		Py = (EN_AND_PR_counter.p_Y/(Lx/4.0)/Ly*1e23/N_a);
		id_P = R*temperature*nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a/1000.0;
		Press = id_P + (Px + Py)/2.0;

		////////////// Trying to change the state ////////////////////////////
		if (RanGen.Random() < 0.5)
			{
				unit_cell_size_variation = true;
				trial_delta_uc = delta_uc * RanGen.Random();
			}
			else
				{
					unit_cell_size_variation = false;
					param_number = RanGen.IRandom(3, params_amount);
					trial_delta_dist_and_angle = (2 * delta_dist_and_angle * RanGen.Random() - delta_dist_and_angle);
				}

		if (unit_cell_size_variation)
			{
				if (Press < 0) // We should decrease the unit cell size
					{
						if (Px < Py) // We should decrease lx
							{
								x_uc = x_uc - trial_delta_uc;
								Lx =  3*x_uc;
							}
							else // We should decrease ly
								{
									y_uc = y_uc - trial_delta_uc;
									Ly = 3*y_uc;
								}
					}
					else // We should increase the unit cell size
						{
							if (Px > Py) // We should increase lx
								{
									x_uc = x_uc - trial_delta_uc;
									Lx =  3*x_uc;
								}
								else // We should increase ly
									{
										y_uc = y_uc - trial_delta_uc;
										Ly =  3*y_uc;
									}
						}

				// Recalculte the coordinates for the expanded or compressed unit cell
				for (int i = 0; i < N; i++)
					{
						coordinates[i].x = coordinates[i].x*(Lx/default_Lx);
						coordinates[i].y = coordinates[i].y*(Ly/default_Ly);
						coordinates[i].x = PBC2D(Lx, coordinates[i].x);
						coordinates[i].y = PBC2D(Ly, coordinates[i].y);
						coordinates[i].damping_coeff = 1.0;
						coordinates[i].ex_field_coeff = empty_field;
						coordinates[i].stat_weight = 1.0;
						coordinates[i].en_and_pr = empty_field;
						// Angles of molecules orientation
						// do not change when expand or compress the unit cell
					}

				//	Realculate the energy of new state
				for(int molA = 0; molA < (N - 1); molA++)
				{
					for(int molB = (molA + 1); molB < N; molB++)
						{
							en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
							en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
							coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
							coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
						}
				}
				weighted_averages_in_central_cell(coordinates, N, Lx, Ly);
				energy_decision = (temp_energy > EN_AND_PR_counter.energy);

				if (energy_decision && !HC_radius && EN_AND_PR_counter.energy < E_INF)
					{

						temp_energy = EN_AND_PR_counter.energy;
						Px = (EN_AND_PR_counter.p_X/(Lx/4.0)/Ly*1e23/N_a);
						Py = (EN_AND_PR_counter.p_Y/(Lx/4.0)/Ly*1e23/N_a);
						id_P = R*temperature*nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a/1000.0;
						Press = id_P + (Px + Py)/2.0;

						// Modification of the unit cell parameters
						if (x_uc - params[1] == 0 && y_uc - params[2] ==0) {cout << "dx: " << x_uc - params[1] << "\t" << "dy: " << y_uc - params[2] << endl;}
						params[1] = x_uc;
						params[2] = y_uc;
						params[3] = sqrt(coordinates[0].x*coordinates[0].x + coordinates[0].y*coordinates[0].y);
						params[4] = acos(coordinates[0].x/params[3])/PI*180.0;
						for (int i = 1; i < params[0]; i++)
							{
								params[3+i*3] = sqrt((coordinates[i].x - coordinates[i-1].x)*(coordinates[i].x - coordinates[i-1].x) + (coordinates[i].y - coordinates[i-1].y)*(coordinates[i].y - coordinates[i-1].y));
								params[4+i*3] = acos((coordinates[i].x - coordinates[i-1].x)/params[3+i*3])/PI*180.0;
							}
						step++;
						cout << step << "\t" << counter  << "\t" << trial_delta_uc << "\t" << trial_delta_dist_and_angle << "\t"
						<< Lx << "\t" << Ly << "\t"
						<< nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a << "\t" << EN_AND_PR_counter.energy / 1000.0 / nPart_in_central_cell << "\t" << Press << "\t" << Px << "\t" << Py <<  endl;
						counter = 0;
					}
					else {counter++;}
			}
			else
				{
					// Shift a molecule in the unit cell by cahnging one of the graph parameters
					params[param_number] += trial_delta_dist_and_angle;
					int molecules = 0;
					//filling the unit cell for the new state
					for(int i = 0; i < 3; i++)
					{
						for(int j = 0; j < 3; j++)
						{
							coordinates[molecules].x = i * x_uc + params[3] * cos(params[4] / 180.0 * PI);
							coordinates[molecules].y = j * y_uc + params[3] * sin(params[4] / 180.0 * PI);
							coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
							coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
							coordinates[molecules].phi = params[5];
							coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
							coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
							coordinates[molecules].damping_coeff = 1.0;
							coordinates[molecules].ex_field_coeff = empty_field;
							coordinates[molecules].stat_weight = 1.0;
							coordinates[molecules].en_and_pr = empty_field;
							molecules++;
							for (int mol = 1; mol < params[0]; mol++)
							{
								int number = 3 + mol * 3;
								coordinates[molecules].x = coordinates[molecules - 1].x + params[number] * cos(params[number + 1] / 180.0 * PI);
								coordinates[molecules].y = coordinates[molecules - 1].y + params[number] * sin(params[number + 1] / 180.0 * PI);
								coordinates[molecules].x = PBC2D(Lx, coordinates[molecules].x);
								coordinates[molecules].y = PBC2D(Ly, coordinates[molecules].y);
								coordinates[molecules].phi = params[number + 2];
								coordinates[molecules].sin_phi = sin(coordinates[molecules].phi / 180.0 * PI);
								coordinates[molecules].cos_phi = cos(coordinates[molecules].phi / 180.0 * PI);
								coordinates[molecules].damping_coeff = 1.0;
								coordinates[molecules].ex_field_coeff = empty_field;
								coordinates[molecules].stat_weight = 1.0;
								coordinates[molecules].en_and_pr = empty_field;
								molecules++;
							}
						}
					}
					// Recalculte the energy of new state
					for(int molA = 0; molA < (N - 1); molA++)
					{
						for(int molB = (molA + 1); molB < N; molB++)
							{
								en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, true);
								en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
								coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
								coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
							}
					}
					weighted_averages_in_central_cell(coordinates, N, Lx, Ly);
					energy_decision = (temp_energy > EN_AND_PR_counter.energy);

					if (energy_decision && !HC_radius && EN_AND_PR_counter.energy < E_INF)
						{
							temp_energy = EN_AND_PR_counter.energy;
							Px = (EN_AND_PR_counter.p_X/(Lx/4.0)/Ly*1e23/N_a);
							Py = (EN_AND_PR_counter.p_Y/(Lx/4.0)/Ly*1e23/N_a);
							id_P = R*temperature*nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a/1000.0;
							Press = id_P + (Px + Py)/2.0;
							step++;
							cout << step << "\t" << counter << "\t" << trial_delta_uc << "\t" << trial_delta_dist_and_angle << "\t"
							<< Lx << "\t" << Ly << "\t"
							<< nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a << "\t" << EN_AND_PR_counter.energy / 1000.0 / nPart_in_central_cell << "\t" << Press << "\t" << Px << "\t" << Py <<  endl;
							counter = 0;
						}
						else {params[param_number] -= trial_delta_dist_and_angle; counter++;}
				}
	}

  cout << "Final params: " << endl;
  for (int i = 0; i < params.size(); i++)
  {
    cout << "Number " << i << ": " << params[i] << endl;
  }

  generate_elongated_cell(params, coordinates, Lx, Ly);
  double E_INF = temp_E_INF;
}
