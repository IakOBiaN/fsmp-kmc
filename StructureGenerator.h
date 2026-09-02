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

  int needed = (int)params[0] * uc_in_x * uc_in_y;
  if (needed > (int)coordinates.size())
  {
    cerr << "ERROR: structure needs " << needed << " molecules, but the coordinates buffer holds only "
         << coordinates.size() << ". Increase its capacity in program_body.cpp (vector<state> coordinates)." << endl;
    exit(1);
  }

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
  }

  // The first unit cell of the final coordinates defines the site lattice of
  // the stabilization mask (no-op unless stabilization_mask = true)
  mask_build(params, coordinates);

  for (int i = 0; i < molecules; i++)
  {
    coordinates[i].damping_coeff = damping_field(coordinates[i].x, Lx); // Lambda^1/2
    coordinates[i].ex_field_coeff = external_field_and_mask(coordinates[i].x, coordinates[i].y, Lx); // u_ext + mask
    coordinates[i].stat_weight = weights_for_central_cell (coordinates[i].x, Lx);
  }

  cout << endl << "Elongated cell was generated: " << endl;
  cout << "N: " << molecules << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
}

void generate_structure(vector <double> &params, string structure_name, vector <state> &coordinates, double &Lx, double &Ly)
{
  if (structure_name == "TMA_HCP_simple_2020")
{
  unit_cell_params.push_back(2);
  unit_cell_params.push_back(11.1);
  unit_cell_params.push_back(19.2258);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(90);
  unit_cell_params.push_back(11.1);
  unit_cell_params.push_back(60.0003);
  unit_cell_params.push_back(90);
}
  // Chicken-wire (honeycomb) cell for the simplified 2020 model: E = -49.56
  // kJ/mol, density 1.064 umol/m2 at T = 0. Metastable in this model, so
  // run it with stabilization_mask = true.
  if (structure_name == "TMA_CW_simple_2020")
{
  unit_cell_params.push_back(4);
  unit_cell_params.push_back(18.9835);
  unit_cell_params.push_back(32.8801);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(29.9997);
  unit_cell_params.push_back(10.9602);
  unit_cell_params.push_back(29.9997);
  unit_cell_params.push_back(89.9995);
  unit_cell_params.push_back(10.96);
  unit_cell_params.push_back(90.0002);
  unit_cell_params.push_back(30.0002);
  unit_cell_params.push_back(10.9602);
  unit_cell_params.push_back(150);
  unit_cell_params.push_back(90.0001);
}
  if (structure_name == "TMA_FL2_simple_2020")
{
  unit_cell_params.push_back(24);
  unit_cell_params.push_back(41.1867);
  unit_cell_params.push_back(71.3478);
  unit_cell_params.push_back(9.28874);
  unit_cell_params.push_back(45);
  unit_cell_params.push_back(89.9474);
  unit_cell_params.push_back(11.1117);
  unit_cell_params.push_back(0.0783955);
  unit_cell_params.push_back(90.0303);
  unit_cell_params.push_back(11.0802);
  unit_cell_params.push_back(120.109);
  unit_cell_params.push_back(90.1096);
  unit_cell_params.push_back(11.1178);
  unit_cell_params.push_back(0.297238);
  unit_cell_params.push_back(90.2973);
  unit_cell_params.push_back(15.6402);
  unit_cell_params.push_back(135.81);
  unit_cell_params.push_back(30.5047);
  unit_cell_params.push_back(11.1392);
  unit_cell_params.push_back(0.392199);
  unit_cell_params.push_back(30.394);
  unit_cell_params.push_back(19.057);
  unit_cell_params.push_back(-59.9215);
  unit_cell_params.push_back(29.811);
  unit_cell_params.push_back(10.9601);
  unit_cell_params.push_back(29.9136);
  unit_cell_params.push_back(90.6418);
  unit_cell_params.push_back(10.9614);
  unit_cell_params.push_back(90.6417);
  unit_cell_params.push_back(29.9372);
  unit_cell_params.push_back(10.9612);
  unit_cell_params.push_back(149.516);
  unit_cell_params.push_back(89.5164);
  unit_cell_params.push_back(15.6136);
  unit_cell_params.push_back(164.98);
  unit_cell_params.push_back(30.0546);
  unit_cell_params.push_back(11.0907);
  unit_cell_params.push_back(179.919);
  unit_cell_params.push_back(29.8062);
  unit_cell_params.push_back(11.1199);
  unit_cell_params.push_back(59.8066);
  unit_cell_params.push_back(30.2144);
  unit_cell_params.push_back(15.5407);
  unit_cell_params.push_back(-15.3172);
  unit_cell_params.push_back(90.1936);
  unit_cell_params.push_back(11.1325);
  unit_cell_params.push_back(0.0313265);
  unit_cell_params.push_back(89.9321);
  unit_cell_params.push_back(19.2179);
  unit_cell_params.push_back(149.95);
  unit_cell_params.push_back(89.7438);
  unit_cell_params.push_back(11.0999);
  unit_cell_params.push_back(-0.129291);
  unit_cell_params.push_back(89.7815);
  unit_cell_params.push_back(11.1);
  unit_cell_params.push_back(-0.0898669);
  unit_cell_params.push_back(89.8569);
  unit_cell_params.push_back(21.9915);
  unit_cell_params.push_back(30.0337);
  unit_cell_params.push_back(29.8465);
  unit_cell_params.push_back(11.0899);
  unit_cell_params.push_back(-0.153647);
  unit_cell_params.push_back(29.8451);
  unit_cell_params.push_back(11.0856);
  unit_cell_params.push_back(-0.0879708);
  unit_cell_params.push_back(29.8566);
  unit_cell_params.push_back(19.2019);
  unit_cell_params.push_back(149.865);
  unit_cell_params.push_back(30.023);
  unit_cell_params.push_back(11.0829);
  unit_cell_params.push_back(-0.0615314);
  unit_cell_params.push_back(29.947);
  unit_cell_params.push_back(15.5846);
  unit_cell_params.push_back(-15.3659);
  unit_cell_params.push_back(90.1694);
}
  //old optimization
	if (structure_name == "TMA_fCW_qB3LYP_PBE_Dreiding_Dhb5.4")
	{
    unit_cell_params.push_back(6);
    unit_cell_params.push_back(17.2328);
    unit_cell_params.push_back(29.8526);
    unit_cell_params.push_back(-1.43533);
    unit_cell_params.push_back(0.294197);
    unit_cell_params.push_back(32.9886);
    unit_cell_params.push_back(9.94999);
    unit_cell_params.push_back(30.0169);
    unit_cell_params.push_back(93.4888);
    unit_cell_params.push_back(9.89061);
    unit_cell_params.push_back(149.3);
    unit_cell_params.push_back(101.053);
    unit_cell_params.push_back(9.88199);
    unit_cell_params.push_back(29.5589);
    unit_cell_params.push_back(32.9726);
    unit_cell_params.push_back(9.94001);
    unit_cell_params.push_back(149.96);
    unit_cell_params.push_back(93.6463);
    unit_cell_params.push_back(9.42501);
    unit_cell_params.push_back(20.3877);
    unit_cell_params.push_back(89.8877);
	}
  //optimized
  if (structure_name == "TMA_CW_qB3LYP_PBE_Dreiding_Dhb5.4")
  {
     unit_cell_params.push_back(4);
     unit_cell_params.push_back(17.2603);
     unit_cell_params.push_back(29.8952);

     unit_cell_params.push_back(0.235947);
     unit_cell_params.push_back(-0.764381);
     unit_cell_params.push_back(26.9998);

     unit_cell_params.push_back(9.9654);
     unit_cell_params.push_back(29.9989);
     unit_cell_params.push_back(86.9999);

     unit_cell_params.push_back(9.96507);
     unit_cell_params.push_back(89.9999);
     unit_cell_params.push_back(27.0004);

     unit_cell_params.push_back(9.96543);
     unit_cell_params.push_back(150.001);
     unit_cell_params.push_back(87.0006);
  }
  if (structure_name == "TMA_CW_aimnet2")
{
  unit_cell_params.push_back(4);
  unit_cell_params.push_back(16.7522);
  unit_cell_params.push_back(29.0156);
  unit_cell_params.push_back(0.361861);
  unit_cell_params.push_back(50.1356);
  unit_cell_params.push_back(33.0001);
  unit_cell_params.push_back(9.67166);
  unit_cell_params.push_back(150);
  unit_cell_params.push_back(93.0012);
  unit_cell_params.push_back(9.67205);
  unit_cell_params.push_back(90.002);
  unit_cell_params.push_back(33.0008);
  unit_cell_params.push_back(9.67171);
  unit_cell_params.push_back(150);
  unit_cell_params.push_back(93.0001);
}
  //optimized
  if (structure_name == "TMA_SF_qB3LYP_PBE_Dreiding_Dhb5.4")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.73934);
     unit_cell_params.push_back(16.869);
     unit_cell_params.push_back(-1.44277);
     unit_cell_params.push_back(1.15426);
     unit_cell_params.push_back(39.4999);
     unit_cell_params.push_back(9.73936);
     unit_cell_params.push_back(60);
     unit_cell_params.push_back(39.5);
  }
  //optimized
  if (structure_name == "TPA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.9);
     unit_cell_params.push_back(11.6897);
     unit_cell_params.push_back(0.770517);
     unit_cell_params.push_back(-0.837901);
     unit_cell_params.push_back(-1.5);
     unit_cell_params.push_back(7.57459);
     unit_cell_params.push_back(50.4947);
     unit_cell_params.push_back(1.00003);
  }
  //optimized
  if (structure_name == "TPA_vertical_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(11.6935);
     unit_cell_params.push_back(9.9);
     unit_cell_params.push_back(-0.740467);
     unit_cell_params.push_back(91.1906);
     unit_cell_params.push_back(88.5061);
     unit_cell_params.push_back(7.5782);
     unit_cell_params.push_back(140.506);
     unit_cell_params.push_back(91.0062);
  }
  if (structure_name == "TPA_horizontal_ladder_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(9.85001);
     unit_cell_params.push_back(14.9538);
     unit_cell_params.push_back(0.459786);
     unit_cell_params.push_back(2.47113);
     unit_cell_params.push_back(3.00003);
     unit_cell_params.push_back(9.35784);
     unit_cell_params.push_back(51.7222);
     unit_cell_params.push_back(79.5001);
  }
  if (structure_name == "TPA_vertical_ladder_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(14.9387);
     unit_cell_params.push_back(9.85001);
     unit_cell_params.push_back(3.1566);
     unit_cell_params.push_back(0.749759);
     unit_cell_params.push_back(93.0001);
     unit_cell_params.push_back(9.54028);
     unit_cell_params.push_back(37.2721);
     unit_cell_params.push_back(-10.4999);
  }
  //optimized
  if (structure_name == "IPA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
  {
     unit_cell_params.push_back(2);
     unit_cell_params.push_back(16.9262);
     unit_cell_params.push_back(7.33665);
     unit_cell_params.push_back(3.16404);
     unit_cell_params.push_back(-2.30584);
     unit_cell_params.push_back(147.971);
     unit_cell_params.push_back(9.85001);
     unit_cell_params.push_back(30.9713);
     unit_cell_params.push_back(328.471);
  }
	//optimized
	if (structure_name == "IPA_vertical_chain_qPBE_Dreiding_Dhb5.0")
	{
		 unit_cell_params.push_back(2);
		 unit_cell_params.push_back(7.49);
		 unit_cell_params.push_back(16.3537);
		 unit_cell_params.push_back(0.472783);
		 unit_cell_params.push_back(-0.646095);
		 unit_cell_params.push_back(58.999);
		 unit_cell_params.push_back(9.80001);
		 unit_cell_params.push_back(56.4236);
		 unit_cell_params.push_back(240.423);
	}
	// optimized
	if (structure_name == "IPA_hexagones_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(24);
		unit_cell_params.push_back(42.6024);
		unit_cell_params.push_back(49.1144);
		//1
		unit_cell_params.push_back(11.4812);
		unit_cell_params.push_back(0.0714106);
		unit_cell_params.push_back(238.318);
		//2
		unit_cell_params.push_back(9.87994);
		unit_cell_params.push_back(120.165);
		unit_cell_params.push_back(297.665);
		//3
		unit_cell_params.push_back(7.59915);
		unit_cell_params.push_back(90.1868);
		unit_cell_params.push_back(178.4);
		//4
		unit_cell_params.push_back(9.88002);
		unit_cell_params.push_back(59.9001);
		unit_cell_params.push_back(237.656);
		//5
		unit_cell_params.push_back(9.87003);
		unit_cell_params.push_back(120.157);
		unit_cell_params.push_back(298.199);
		//6
		unit_cell_params.push_back(7.34052);
		unit_cell_params.push_back(89.7872);
		unit_cell_params.push_back(177.856);
		//7
		unit_cell_params.push_back(12.4172);
		unit_cell_params.push_back(22.7796);
		unit_cell_params.push_back(-1.50074);
		//8
		unit_cell_params.push_back(9.88001);
		unit_cell_params.push_back(240.147);
		unit_cell_params.push_back(58.102);
		//9
		unit_cell_params.push_back(9.87999);
		unit_cell_params.push_back(300.102);
		unit_cell_params.push_back(118.585);
		//10
		unit_cell_params.push_back(7.54303);
		unit_cell_params.push_back(270.085);
		unit_cell_params.push_back(-1.72634);
		//11
		unit_cell_params.push_back(9.87142);
		unit_cell_params.push_back(240.131);
		unit_cell_params.push_back(57.4818);
		//12
		unit_cell_params.push_back(9.87005);
		unit_cell_params.push_back(299.981);
		unit_cell_params.push_back(118.103);
		//13
		unit_cell_params.push_back(9.86024);
		unit_cell_params.push_back(0.777899);
		unit_cell_params.push_back(178.55);
		//14
		unit_cell_params.push_back(9.88015);
		unit_cell_params.push_back(60.0499);
		unit_cell_params.push_back(237.807);
		//15
		unit_cell_params.push_back(9.87447);
		unit_cell_params.push_back(120.041);
		unit_cell_params.push_back(298.541);
		//16
		unit_cell_params.push_back(7.55718);
		unit_cell_params.push_back(90.0471);
		unit_cell_params.push_back(178.355);
		//17
		unit_cell_params.push_back(9.87826);
		unit_cell_params.push_back(60.0484);
		unit_cell_params.push_back(237.732);
		//18
		unit_cell_params.push_back(9.87009);
		unit_cell_params.push_back(120.232);
		unit_cell_params.push_back(298.183);
		//19
		unit_cell_params.push_back(7.48786);
		unit_cell_params.push_back(29.3029);
		unit_cell_params.push_back(57.8625);
		//20
		unit_cell_params.push_back(9.86602);
		unit_cell_params.push_back(299.677);
		unit_cell_params.push_back(117.237);
		//21
		unit_cell_params.push_back(7.47892);
		unit_cell_params.push_back(269.811);
		unit_cell_params.push_back(-1.65904);
		//22
		unit_cell_params.push_back(9.88013);
		unit_cell_params.push_back(239.841);
		unit_cell_params.push_back(57.6999);
		//23
		unit_cell_params.push_back(9.87092);
		unit_cell_params.push_back(300.2);
		unit_cell_params.push_back(118.231);
		//24
		unit_cell_params.push_back(7.41976);
		unit_cell_params.push_back(270.111);
		unit_cell_params.push_back(-1.8894);
	}
	//optimized
	if (structure_name == "PA_horizontal_chain_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(2);
	  unit_cell_params.push_back(12.5557);
	  unit_cell_params.push_back(9.67631);
		unit_cell_params.push_back(0.547972);
		unit_cell_params.push_back(0.0439552);
		unit_cell_params.push_back(127.331);
		unit_cell_params.push_back(9.91771);
		unit_cell_params.push_back(51.2729);
		unit_cell_params.push_back(307.331);
	}
	//optimized
	if (structure_name == "PA_vertical_chain_qPBE_Dreiding_Dhb5.0")
	{
		unit_cell_params.push_back(2);
		unit_cell_params.push_back(9.16996);
		unit_cell_params.push_back(12.9454);
		unit_cell_params.push_back(2.95443);
		unit_cell_params.push_back(0.146893);
		unit_cell_params.push_back(34.7622);
		unit_cell_params.push_back(9.89);
		unit_cell_params.push_back(40.7623);
		unit_cell_params.push_back(214.762);
	}

  generate_elongated_cell(params, coordinates, Lx, Ly);
}

// Build the 3x3 tiling of the unit cell described by params into coordinates and
// return its total energy in J/mol (summed over the 9*params[0] molecules, no
// external fields). Sets the global HC_radius flag if any pair overlaps the hard core.
double optimizer_tiling_energy(vector <double> &params, vector <state> &coordinates, double &Lx, double &Ly, double beta)
{
  results empty_field;
  results en_and_press;
  int N = (int)params[0] * 3 * 3;
  Lx = params[1] * 3;
  Ly = params[2] * 3;
  double x_uc = params[1];
  double y_uc = params[2];
  int molecules = 0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
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

  HC_radius = false;
  for (int molA = 0; molA < (N - 1); molA++)
  {
    for (int molB = (molA + 1); molB < N; molB++)
    {
      en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, false);
      en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
      coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
      coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
    }
  }
  weighted_averages_in_central_cell(coordinates, N, Lx, Ly);
  return EN_AND_PR_counter.energy;
}

// Refine a rough unit cell (structure_name = "calculate"): adaptive random descent
// over the cell sides and the molecular placement parameters. One random parameter
// is perturbed at a time; a move is accepted only when the energy of the 3x3 tiling
// strictly decreases and no hard-core overlap appears. Lengths and angles carry
// separate step sizes; after stall_limit consecutive rejections both steps are
// halved, and convergence is declared once they fall below the thresholds.
//
// params[3] and params[4] (the offset of the first molecule from the cell origin)
// only translate the lattice as a whole, so they are never perturbed.
//
// Deterministic when compiled with -DFSMP_RANDOM_SEED=<n>.
// Stage 0 helper: rescale every length parameter (cell sides, intermolecular
// distances) by a common factor s relative to the starting values and return
// the energy of the rescaled cell.
static double scaled_cell_energy(vector <double> &params, const vector <double> &start,
                                 const vector<int> &dof_scale, double s,
                                 vector <state> &coordinates, double &Lx, double &Ly, double beta)
{
  for (size_t i = 0; i < dof_scale.size(); i++)
  {
    params[dof_scale[i]] = start[dof_scale[i]] * s;
  }
  return optimizer_tiling_energy(params, coordinates, Lx, Ly, beta);
}

void generate_structure(vector <double> &params, vector <state> &coordinates, double &Lx, double &Ly)
{
  double temp_E_INF = E_INF;
  E_INF = 1e200;   // expose the true repulsion to the optimizer instead of the flat cap
  double beta = 1.0 / (R * 300);
  int n_params = (int)params.size();
  int N = (int)params[0] * 3 * 3;

  // Optimizable degrees of freedom, split by type: lengths in A, angles in degrees
  vector<int> dof_len, dof_ang;
  dof_len.push_back(1);   // cell side x
  dof_len.push_back(2);   // cell side y
  if (n_params > 5) { dof_ang.push_back(5); }   // rotation of the first molecule
  for (int mol = 1; 3 + mol * 3 + 2 < n_params; mol++)
  {
    dof_len.push_back(3 + mol * 3);       // distance from the previous molecule
    dof_ang.push_back(3 + mol * 3 + 1);   // direction of that displacement
    dof_ang.push_back(3 + mol * 3 + 2);   // own rotation
  }
  int n_dof = (int)(dof_len.size() + dof_ang.size());

  double step_len = 0.3;               // A
  double step_ang = 5.0;               // deg
  const double step_len_min = 0.001;   // convergence thresholds
  const double step_ang_min = 0.02;
  const int  stall_limit = 300;        // consecutive rejections before halving the steps
  const long iter_cap = 2000000;       // safety net

  cout << endl << "Unit cell optimization: " << n_dof << " degrees of freedom" << endl;

  // --- Stage 0: uniform scaling of the whole cell ---------------------------
  // One common factor rescales the cell sides and all intermolecular distances
  // (the first-molecule offset too: it only shifts the lattice). The energy
  // along this ray is scanned globally from the hard-core edge outwards: a
  // local search is blind on the two flat plateaus of the landscape (the
  // capped repulsion around the core and the exact zero beyond the cutoff),
  // and the scan is cheap. The deepest bound scale wins and is refined by a
  // deterministic 1D pattern search. When the whole ray is repulsive (bound
  // states need different molecule orientations - rotations come with stage
  // 1), the cell starts dense instead, so the stage-1 descent has angular
  // gradients to work with rather than the flat zero landscape beyond the
  // cutoff, from which no local move is ever accepted.
  vector<int> dof_scale = dof_len;
  dof_scale.push_back(3);
  vector<double> start = params;
  double s = 1.0;
  double energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);

  // the smallest overlap-free scale (pair distances grow linearly with s, so
  // the hard-core edge is monotone): grow out of an overlap, or walk down to
  // the edge from a loose start
  double s_lo = 1.0;
  if (HC_radius)
  {
    cout << "The starting cell has hard-core overlaps; growing it" << endl;
    bool separated = false;
    for (int i = 0; i < 64 && !separated; i++)
    {
      s_lo *= 1.1;
      scaled_cell_energy(params, start, dof_scale, s_lo, coordinates, Lx, Ly, beta);
      separated = !HC_radius;
    }
    if (!separated)
    {
      cerr << "ERROR: molecules still overlap after growing the cell hundreds-fold; "
           << "the unit_cell geometry is degenerate (coinciding molecules?)." << endl;
      exit(1);
    }
  }
  else
  {
    while (s_lo > 0.02)
    {
      scaled_cell_energy(params, start, dof_scale, s_lo * 0.9, coordinates, Lx, Ly, beta);
      if (HC_radius) { break; }
      s_lo *= 0.9;
    }
  }

  // geometric scan of the ray: the deepest bound scale, plus the densest
  // sane scale (clear of the capped wall) as the fallback start
  const double dense_bar = 30000.0;   // J/mol per molecule
  double best_s = -1.0, best_e = 0.0, dense_s = -1.0;
  s = s_lo;
  for (int i = 0; i < 120; i++)
  {
    energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);
    if (!HC_radius)
    {
      if (best_s < 0 || energy < best_e) { best_s = s; best_e = energy; }
      if (dense_s < 0 && energy / nPart_in_central_cell < dense_bar) { dense_s = s; }
      if (energy == 0.0) { break; }   // every pair is beyond the cutoff
    }
    s *= 1.05;
  }
  if (best_e < 0.0)
  {
    s = best_s;
  }
  else if (dense_s > 0)
  {
    s = dense_s;
    cout << "No bound cell size along the scaling ray (the starting orientations "
         << "repel at every distance); starting dense for the angular search." << endl;
  }
  else
  {
    cerr << "ERROR: no workable cell size: the cell is strongly repulsive at "
         << "every scale (coinciding molecules?)." << endl;
    exit(1);
  }
  energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);

  // local refinement of the scale; pointless on a repulsive ray, where it
  // would only dilute the cell back towards the zero plateau
  if (best_e < 0.0)
  {
    double h = 0.05;                   // step of the pattern search on the scale factor
    long scale_evals = 0;
    while (h > 0.001 && scale_evals < 10000)
    {
      double e_up = scaled_cell_energy(params, start, dof_scale, s + h, coordinates, Lx, Ly, beta);
      bool up_ok = !HC_radius && e_up < energy;
      double e_dn = 0;
      bool dn_ok = false;
      if (s - h > 0.01)
      {
        e_dn = scaled_cell_energy(params, start, dof_scale, s - h, coordinates, Lx, Ly, beta);
        dn_ok = !HC_radius && e_dn < energy;
      }
      scale_evals += 2;
      if (up_ok && (!dn_ok || e_up <= e_dn)) { s += h; energy = e_up; }
      else if (dn_ok)                        { s -= h; energy = e_dn; }
      else                                   { h /= 2.0; }
    }
    energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);
  }
  HC_radius = false;
  cout << "Cell scaling: factor " << s << ", energy " << energy / 1000.0 / nPart_in_central_cell
       << " kJ/mol per molecule" << endl;

  // --- Stage 1: adaptive random descent over the individual parameters ------
  cout << "Density: " << nPart_in_central_cell * (1.0e+26) / (Lx * Ly) / N_a << "\t"
       << " Energy: " << energy / 1000.0 / nPart_in_central_cell << endl;
  write_xyz_file(unit_cell_name, N, density, Lx, Ly, temperature, coordinates, 0, 1, true);

  long iter = 0, accepted = 0;
  int stall = 0;
  while (iter < iter_cap)
  {
    iter++;
    int pick = RanGen.IRandom(0, n_dof - 1);
    int param_number;
    double step;
    if (pick < (int)dof_len.size()) { param_number = dof_len[pick]; step = step_len; }
    else { param_number = dof_ang[pick - dof_len.size()]; step = step_ang; }

    double old_value = params[param_number];
    params[param_number] = old_value + (2.0 * RanGen.Random() - 1.0) * step;

    bool bad = (step == step_len && params[param_number] <= 0.0);  // lengths must stay positive
    double trial = bad ? 0.0 : optimizer_tiling_energy(params, coordinates, Lx, Ly, beta);
    if (!bad && trial < energy && !HC_radius)
    {
      energy = trial;
      accepted++;
      stall = 0;
      write_xyz_file(unit_cell_name, N, density, Lx, Ly, temperature, coordinates, 0, 1, false);
      cout << "Density: " << nPart_in_central_cell * (1.0e+26) / (Lx * Ly) / N_a << "\t"
           << " Energy: " << energy / 1000.0 / nPart_in_central_cell << endl;
    }
    else
    {
      params[param_number] = old_value;
      stall++;
      if (stall >= stall_limit)
      {
        step_len /= 2.0;
        step_ang /= 2.0;
        stall = 0;
        if (step_len < step_len_min && step_ang < step_ang_min) { break; }
        cout << "Steps halved to " << step_len << " A / " << step_ang
             << " deg (iteration " << iter << ", accepted " << accepted << ")" << endl;
      }
    }
  }

  // Rebuild at the accepted parameters and report
  energy = optimizer_tiling_energy(params, coordinates, Lx, Ly, beta);
  HC_radius = false;
  cout << "Optimization " << ((iter < iter_cap) ? "converged" : "hit the iteration cap")
       << " after " << iter << " iterations (" << accepted << " accepted)" << endl;
  cout << "Final energy per molecule: " << energy / 1000.0 / nPart_in_central_cell << " kJ/mol" << endl;
  cout << "Final density: " << nPart_in_central_cell * (1.0e+26) / (Lx * Ly) / N_a << " mkmol/m2" << endl;
  cout << "Final params: " << endl;
  streamsize params_prec = cout.precision(17);
  for (size_t i = 0; i < params.size(); i++)
  {
    cout << "Number " << i << ": " << params[i] << endl;
  }
  cout.precision(params_prec);

  generate_elongated_cell(params, coordinates, Lx, Ly);
  E_INF = temp_E_INF;   // restore the hard-core cap disabled for the optimization
}
