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
    coordinates[i].damping_coeff = damping_field(coordinates[i].x, Lx); // Lambda^1/2
    coordinates[i].ex_field_coeff = external_field(coordinates[i].x, Lx); // u_ext
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
  unit_cell_params.push_back(11.100);
  unit_cell_params.push_back(19.200);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(0);
  unit_cell_params.push_back(90);
  unit_cell_params.push_back(11.089); 
  unit_cell_params.push_back(59.967); 
  unit_cell_params.push_back(90);
}
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
  // (the first-molecule offset too: it only shifts the lattice). A start with
  // hard-core overlaps is first grown until the overlaps disappear, then a 1D
  // pattern search (deterministic, no randomness) finds the optimal overall
  // size. Only after that the per-parameter refinement below takes over.
  vector<int> dof_scale = dof_len;
  dof_scale.push_back(3);
  vector<double> start = params;
  double s = 1.0;
  double energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);
  if (HC_radius)
  {
    cout << "The starting cell has hard-core overlaps; growing it" << endl;
    bool separated = false;
    for (int i = 0; i < 64 && !separated; i++)
    {
      s *= 1.1;
      energy = scaled_cell_energy(params, start, dof_scale, s, coordinates, Lx, Ly, beta);
      separated = !HC_radius;
    }
    if (!separated)
    {
      cerr << "ERROR: molecules still overlap after growing the cell hundreds-fold; "
           << "the unit_cell geometry is degenerate (coinciding molecules?)." << endl;
      exit(1);
    }
  }
  double h = 0.05;                     // step of the pattern search on the scale factor
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
  for (size_t i = 0; i < params.size(); i++)
  {
    cout << "Number " << i << ": " << params[i] << endl;
  }

  generate_elongated_cell(params, coordinates, Lx, Ly);
  E_INF = temp_E_INF;   // restore the hard-core cap disabled for the optimization
}
