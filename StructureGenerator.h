using namespace std;

void generate_structure(string structure_name, int &nPart, double &density, vector <state> &coordinates, double &Lx, double &Ly, double state_dens)
{
	if (structure_name == "fCW")
	{
		int cells = nPart/6.0/8.0;  //4 molecules in unit cell (amount)
		int number_in_x = sqrt(cells) + 1;  //unit cells along x directions
		int number_in_y = number_in_x;  //unit cells along y directions
		nPart = number_in_x*number_in_y*48.0;
		double state_Ly = sqrt((1.0e+26)*nPart/8.0/(state_dens*N_a/sqrt(3.0)));  //y-size of the cell from density
		double y_uc = state_Ly/number_in_y;  //y-size of the unit cell
		double x_uc = y_uc/sqrt(3.0);  //x-size of the unit cell
		double h_bond_dist = y_uc/3.0;  //hydrogen bond distance related
		number_in_x *= 8;

		Lx = number_in_x*x_uc; //add vacuum slab at both sides
		Ly = number_in_y*y_uc;  //y-size of the cell

		int molecule = 0; // Molecules counter
		for(int i = 0; i < number_in_x; i++)
			{
			 for(int j = 0; j < number_in_y; j++)
				{
				  //unit cell
				 coordinates[molecule].x = i*x_uc;
				  coordinates[molecule].y = j*y_uc;
				  coordinates[molecule].phi = 30.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;

				  coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
				  coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
				  coordinates[molecule].phi = 90.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;

				  coordinates[molecule].x = coordinates[molecule-1].x;
				  coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
				  coordinates[molecule].phi = 30.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;

				  coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
				  coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
				  coordinates[molecule].phi = 90.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;

				  coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
				  coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist / 2.0;
				  coordinates[molecule].phi = 90.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;

				  coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
				  coordinates[molecule].y = coordinates[molecule-1].y - h_bond_dist*sin(30.0/180.0*PI) - h_bond_dist;
				  coordinates[molecule].phi = 90.0;
				  coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
				  coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
							coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
							coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
							coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
				  molecule++;
			  }
			}

			for (int i = 0; i < molecule; i++)
				{
					double abs_x;
					if (coordinates[i].x > Lx/2.0){abs_x = coordinates[i].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[i].x;}
					double ksi = 32.0*abs_x/Lx;
					if (ksi > 8.0)
					{
						coordinates[i] = coordinates[molecule-1];
						molecule--;
						i--;
					}
				}
				nPart = molecule;

		 cout << endl << "Filled Chichen-Wire TMA Structure in central cell: " << endl;
		 cout << "N: " << molecule << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
	}
//////////////////////////////////////////////////////////////////////////////////
  if (structure_name == "CW")
  {
    int cells = nPart/4.0/8.0;  //4 molecules in unit cell (amount)
    int number_in_x = sqrt(cells) + 1;  //unit cells along x directions
    int number_in_y = number_in_x;  //unit cells along y directions
    nPart = number_in_x*number_in_y*32.0;
    double state_Ly = sqrt((1.0e+26)*nPart/8.0/(state_dens*N_a/sqrt(3.0)));  //y-size of the cell from density
    double y_uc = state_Ly/number_in_y;  //y-size of the unit cell
    double x_uc = y_uc/sqrt(3.0);  //x-size of the unit cell
    double h_bond_dist = y_uc/3.0;  //hydrogen bond distance related
    number_in_x *= 8;
    cout << "AAA:" << x_uc << " " << y_uc << " " << h_bond_dist << endl;

    Lx = number_in_x*x_uc; //add vacuum slab at both sides
    Ly = number_in_y*y_uc;  //y-size of the cell

    int molecule = 0; // Molecules counter
    for(int i = 0; i < number_in_x; i++)
        {
         for(int j = 0; j < number_in_y; j++)
            {
              //unit cell
             coordinates[molecule].x = i*x_uc;
              coordinates[molecule].y = j*y_uc;
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
              molecule++;

              coordinates[molecule].x = coordinates[molecule-1].x + h_bond_dist*cos(30.0/180.0*PI);
              coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
              coordinates[molecule].phi = 90.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
              molecule++;

              coordinates[molecule].x = coordinates[molecule-1].x;
              coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist;
              coordinates[molecule].phi = 30.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
              molecule++;

              coordinates[molecule].x = coordinates[molecule-1].x - h_bond_dist*cos(30.0/180.0*PI);
              coordinates[molecule].y = coordinates[molecule-1].y + h_bond_dist*sin(30.0/180.0*PI);
              coordinates[molecule].phi = 90.0;
              coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
              coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    					coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    					coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    					coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
              molecule++;
          }
        }

        for (int i = 0; i < molecule; i++)
    		{
    			double abs_x;
    			if (coordinates[i].x > Lx/2.0){abs_x = coordinates[i].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[i].x;}
    			double ksi = 32.0*abs_x/Lx;
    			if (ksi > 8.0)
    			{
    				coordinates[i] = coordinates[molecule-1];
    				molecule--;
    				i--;
    			}
    		}
    		nPart = molecule;

     cout << endl << "Chichen-Wire TMA Structure in central cell: " << endl;
     cout << "N: " << molecule << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
  }
  if (structure_name == "SF")
  {
    int cells = (nPart/2.0/8.0);  //2 molecules in unit cell (amount)
    int number_in_x = sqrt(cells) + 1;  //unit cells along x-directions
    int number_in_y = number_in_x;  //unit cells along y-directions
    nPart = number_in_x*number_in_y*16.0;
    double state_Ly = sqrt((1.0e+26)*nPart/8.0/(state_dens*N_a/sqrt(3)));
    double y_uc = state_Ly/number_in_y;  //y-size of the unit cell
    double x_uc = y_uc/sqrt(3.0);  //x-size of the unit cell
    double h3_bond_dist = x_uc;  //hydrogen 3-bond distance related
    number_in_x *= 8;

    Lx = number_in_x*x_uc;  //initial x-size of the crystal
    Ly = number_in_y*y_uc;  //y-size of the cell
    cout << "Lx: " << Lx << "\t" << "Ly: " << Ly << endl;

    int molecule = 0; // Molecules counter
    for(int i = 0; i < number_in_x; i++)
    		{
         for(int j = 0; j < number_in_y; j++)
    			{
    				coordinates[molecule].x = h3_bond_dist*sin(60.0/180.0*PI)/2.0 + i*x_uc;
    				coordinates[molecule].y = h3_bond_dist/4.0 + j*y_uc;
    				coordinates[molecule].phi = 30.0;
    				coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
    				coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    				coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    				coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    				coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
    				molecule++;

    				coordinates[molecule].x = coordinates[molecule-1].x + h3_bond_dist*cos(60.0/180.0*PI);
    				coordinates[molecule].y = coordinates[molecule-1].y + h3_bond_dist*sin(60.0/180.0*PI);
    				coordinates[molecule].phi = 46.0;
    				coordinates[molecule].sin_phi = sin(coordinates[molecule].phi/180.0*PI);
    				coordinates[molecule].cos_phi = cos(coordinates[molecule].phi/180.0*PI);
    				coordinates[molecule].damping_coeff = damping_field(coordinates[molecule].x, Lx); // Lambda^1/2
    				coordinates[molecule].ex_field_coeff = external_field(coordinates[molecule].x, Lx); // u_ext
    				coordinates[molecule].stat_weight = weights_for_central_cell (coordinates[molecule].x, Lx);
    				molecule++;
    			}
    		}

    		for (int i = 0; i < molecule; i++)
    		{
    			double abs_x;
    			if (coordinates[i].x > Lx/2.0){abs_x = coordinates[i].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[i].x;}
    			double ksi = 32.0*abs_x/Lx;
    			if (ksi > 8.0)
    			{
    				coordinates[i] = coordinates[molecule-1];
    				molecule--;
    				i--;
    			}
    		}
    		nPart = molecule;
    		cout << "nPart after cutting: " << nPart << endl;

    		cout << endl << "SuperFlower TMA Structure in central cell: " << endl;
    	  cout << "N: " << molecule << "\t" << "Lx and Ly in A: " << Lx << " and " << Ly << endl;
  }

}

void generate_structure(vector <double> &params, vector <state> &coordinates, double &Lx, double &Ly)
{
  results empty_field;
  results en_and_press;
  int params_amount = params.size() - 1;
  double delta_uc = 0.5;
  double temp_energy = 1e10;
  double convergency = 1e10;
  bool first = true;

  //while (convergency > 5000)
  for (int aaa = 0; aaa < 1000; aaa++)
  {
    int param_number;
    double temp_delta;
    if (!first)
    {
      param_number = RanGen.IRandom(1, params_amount);
      temp_delta = (2 * delta_uc * RanGen.Random() - delta_uc);
    }
    else
    {
      temp_delta = 0;
    }

    params[param_number] += temp_delta;

    Lx = params[1];
    Ly = params[2];
    coordinates[0].x = params[3] * cos(params[4] / 180.0 * PI);
    coordinates[0].y = params[3] * sin(params[4] / 180.0 * PI);
    coordinates[0].x = PBC2D(Lx, coordinates[0].x);
    coordinates[0].y = PBC2D(Ly, coordinates[0].y);
    coordinates[0].phi = params[5];
    coordinates[0].sin_phi = sin(coordinates[0].phi / 180.0 * PI);
    coordinates[0].cos_phi = cos(coordinates[0].phi / 180.0 * PI);
    coordinates[0].damping_coeff = 1.0;
    coordinates[0].ex_field_coeff = empty_field;
    coordinates[0].stat_weight = 1.0;
    coordinates[0].en_and_pr = empty_field;
    for (int mol = 1; mol < params[0]; mol++)
    {
      int number = 3 + mol * 3;
      coordinates[mol].x = coordinates[mol - 1].x + params[number] * cos(params[number + 1] / 180.0 * PI);
      coordinates[mol].y = coordinates[mol - 1].y + params[number] * sin(params[number + 1] / 180.0 * PI);
      coordinates[mol].x = PBC2D(Lx, coordinates[mol].x);
      coordinates[mol].y = PBC2D(Ly, coordinates[mol].y);
      coordinates[mol].phi = params[number + 2];
      coordinates[mol].sin_phi = sin(coordinates[mol].phi / 180.0 * PI);
      coordinates[mol].cos_phi = cos(coordinates[mol].phi / 180.0 * PI);
      coordinates[mol].damping_coeff = 1.0;
      coordinates[mol].ex_field_coeff = empty_field;
      coordinates[mol].stat_weight = 1.0;
      coordinates[mol].en_and_pr = empty_field;
    }
    double beta = 1.0 / (R * 300);
    int N = params[0];

    for(int molA = 0; molA < (N - 1); molA++)
  	{
  		for(int molB = (molA + 1); molB < N; molB++)
  			{
  				en_and_press = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta, false);
  				en_and_press = en_and_press / 2.0;  //for molecules pair to value per molecule
  				coordinates[molA].en_and_pr = coordinates[molA].en_and_pr + en_and_press;
  				coordinates[molB].en_and_pr = coordinates[molB].en_and_pr + en_and_press;
  			}
    }

    weighted_averages_in_central_cell(coordinates, N, Lx, Ly);
    if (first && HC_radius)
    {
      cout << "AHTUNG!!! HC_radius!!!" << endl;
    }

    if ((temp_energy >= EN_AND_PR_counter.energy) && !HC_radius)
    {
      temp_energy = EN_AND_PR_counter.energy;
      write_xyz_file_TMA (N, density, Lx, Ly, temperature, coordinates, 0, 1, false);
      cout << "YES! Density: " << nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a << "\t" << " Energy: " << EN_AND_PR_counter.energy / 1000.0 / nPart_in_central_cell << endl;
    }
    else
    {
      params[param_number] -= temp_delta;
    }

    cout << "Density: " << nPart_in_central_cell * (1.0e+26) / (Lx*Ly) / N_a << "\t" << " Energy: " << EN_AND_PR_counter.energy / 1000.0 / nPart_in_central_cell << endl;
    first = false;
    HC_radius = false;
  }
}
