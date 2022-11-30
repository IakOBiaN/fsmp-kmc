void centering_the_crystal(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{
	// Calculate x-coordinate of mass center of the system
	double xc = 0;
	int mol_in_right_gas = 0;
	int mol_in_left_gas = 0;
	vector <int> right_end (nPart*4);
	vector <int> left_end (nPart*4);
	for (int i = 0; i < nPart; i++)
		{
			xc += coordinates[i].x;
			if (coordinates[trialPart].x > 7.0*Lx/8.0){right_end[mol_in_right_gas] = i; mol_in_right_gas += 1;}
			if (coordinates[trialPart].x < Lx/8.0){left_end[mol_in_left_gas] = i; mol_in_left_gas += 1;}
		}
	xc = xc/nPart;

	results delta_EP;
	results old_EP;
	results new_EP;
	results old_EP_Part;
	results new_EP_Part;

	bool find = false;
	state new_coordinates;

	if (mol_in_left_gas > 0 && mol_in_right_gas > 0)
	{
		if (xc > Lx/2.0)
			{
				int trialPart = RanGen.IRandom(0,(mol_in_right_gas-1));
				trialPart = right_end[trialPart];
				new_coordinates = coordinates[trialPart];
				new_coordinates.x = RanGen.Random()*Lx/8.0;
				new_coordinates.y = Ly * RanGen.Random();
				charges_coordinates (new_coordinates);
			}
		else
			{
				int trialPart = RanGen.IRandom(0,(mol_in_left_gas-1));
				trialPart = left_end[trialPart];
				new_coordinates = coordinates[trialPart];
				new_coordinates.x = (7.0 + RanGen.Random())*Lx/8.0;
				new_coordinates.y = Ly * RanGen.Random();
				charges_coordinates (new_coordinates);
			}

		for (int l = 0; l < nPart; l++)
			{
				if (l == trialPart){continue;}
				old_EP = old_EP +  energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta, false);
				new_EP = new_EP + energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, false);
			}

		old_EP.energy += coordinates[trialPart].ex_field_coeff.energy;
		new_EP.energy += new_coordinates.ex_field_coeff.energy;
		delta_EP.energy = new_EP.energy - old_EP.energy;

		if(RanGen.Random() < exp(-delta_EP.energy*beta) && !HC_radius)
			{
				for (int l = 0; l < nPart; l++)
					{
						if (l == trialPart){continue;}
						// Here we should calculate the pressure if the change is accepted
						old_EP_Part = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta, true);
						new_EP_Part = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, true);
						delta_EP = (new_EP_Part - old_EP_Part)/2.0;
						coordinates[l].en_and_pr = coordinates[l].en_and_pr + delta_EP;
						new_coordinates.en_and_pr = new_coordinates.en_and_pr + delta_EP;
					}
					new_coordinates.en_and_pr = new_coordinates.en_and_pr + (new_coordinates.ex_field_coeff - coordinates[trialPart].ex_field_coeff);
					coordinates[trialPart] = new_coordinates;     // Update position
			}

	}
}
