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

  string cell_out = unit_cell_name;
  size_t dot = cell_out.find_last_of('.');
  size_t slash = cell_out.find_last_of("/\\");
  if (dot != string::npos && (slash == string::npos || dot > slash)) { cell_out.erase(dot); }
  cell_out += ".cell";
  write_cell_file(cell_out, params, "optimized on " + p_name);
  cout << "Optimized cell written to " << cell_out << endl;

  generate_elongated_cell(params, coordinates, Lx, Ly);
  E_INF = temp_E_INF;   // restore the hard-core cap disabled for the optimization
}
