void Rosenbluth_iteration(double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &dt, double &beta)
{
	int trialPart = 0;

	vector <double> mobility_histogram(nPart);
	vector <double> ln_mob(nPart);
	double s = 0;
	dt = 0;

	// Calculating the log mobilities of all the molecules
	for(int i = 0; i < nPart; i++)
		{
			ln_mob[i] = (2.0*coordinates[i].en_and_pr.energy - coordinates[i].ex_field_coeff.energy)*beta;
//			if (ln_mob[i] > 100){ln_mob[i] = 100;}
			s += ln_mob[i];
		}
	s = s/nPart; // Average mobility of the particles

	// Construct the mobility histogram
	mobility_histogram[0] = exp(ln_mob[0]- s);
	for(int i = 1; i < nPart; i++)
		{
			mobility_histogram[i] = mobility_histogram[i-1] + exp(ln_mob[i] - s);
		}

	// Construct the normalized mobility histogram
	for(int i = 0; i < nPart; i++)
		{
			mobility_histogram[i] = mobility_histogram[i]/mobility_histogram[nPart-1];
		}

	// Random selection of the molecule according to the Rosenbluth scheme
	double Rp = mobility_histogram[nPart-1]*RanGen.Random();
	if(Rp < mobility_histogram[0]) {trialPart = 0;}
		else{
					for(int i = 1; i < nPart; i++)
						{
							if(Rp >= mobility_histogram[i-1] && Rp < mobility_histogram[i]) {trialPart = i; break;}
						}
				}
	state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

	results delta_EP;
	results old_EP;
	results new_EP;
	results old_EP_Part;
	results new_EP_Part;

	// Move and Rotate a molecule
	new_coordinates.x = Lx * RanGen.Random();
	new_coordinates.y = Ly * RanGen.Random();
	new_coordinates.phi = 360.0*RanGen.Random();
	new_coordinates.sin_phi = sin(new_coordinates.phi/180.0*PI);
	new_coordinates.cos_phi = cos(new_coordinates.phi/180.0*PI);
	new_coordinates.damping_coeff = damping_field(new_coordinates.x, Lx); // Lambda^1/2
	new_coordinates.ex_field_coeff = external_field(new_coordinates.x, Lx); // u_ext
	new_coordinates.stat_weight = weights_for_central_cell(new_coordinates.x, Lx);

	// Change the state of the system if molecules do not overlap
	for (int l = 0; l < nPart; l++)
		{
			if (l == trialPart){continue;}
			energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, true);
			if (HC_radius == true) {dt = 0; break;} // If the molecules overlap, the duration of new stste is supposed to be zero
		}

	if (HC_radius == false)
	{
		for (int l = 0; l < nPart; l++)
			{
				if (l == trialPart){continue;}
				//Choose exact or numerical energy and pressure calculation
				old_EP_Part = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly, beta, true);
				new_EP_Part = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, true);
				delta_EP = (new_EP_Part - old_EP_Part)/2.0;
				coordinates[l].en_and_pr = coordinates[l].en_and_pr + delta_EP;
				new_coordinates.en_and_pr = new_coordinates.en_and_pr + delta_EP;
			}
			new_coordinates.en_and_pr = new_coordinates.en_and_pr + (new_coordinates.ex_field_coeff - coordinates[trialPart].ex_field_coeff);
			coordinates[trialPart] = new_coordinates;     // Update position of the trialPart

			// Calculating the log mobilities of all the molecules in new configuration of the system
			s = 0;
			for(int i = 0; i < nPart; i++)
				{
					ln_mob[i] = (2.0*coordinates[i].en_and_pr.energy - coordinates[i].ex_field_coeff.energy)*beta;
//					if (ln_mob[i] > 100){ln_mob[i] = 100;}
					s += ln_mob[i];
				}
			s = s/nPart; // Average mobility of the particles in new configuration of the system
			// Construct the mobility histogram in new configuration of the system
			mobility_histogram[0] = exp(ln_mob[0]- s);
			for(int i = 1; i < nPart; i++)
				{
					mobility_histogram[i] = mobility_histogram[i-1] + exp(ln_mob[i] - s);
				}
			// Calculate the duration of new configuration of the system
			dt = exp(-s)/mobility_histogram[nPart-1];
	}
		HC_radius = false;
}
