double Rosenbluth_iteration(double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &dt, double &beta, int &iter, int &trialPart, bool findTrialPart)
{

	//int trialPart = 0;
	vector <double> mob_histogram(nPart);
	double Rp;
	if (findTrialPart)
	{
		double s = 0;
		for(int mol = 0; mol < nPart; mol++)
		{
			s += (2.0*coordinates[mol].en_and_pr.energy - coordinates[mol].ex_field_coeff.energy)*beta;
		}
		s = s/nPart;

		double sq = 0;
		for(int mol = 0; mol < nPart; mol++)
		{
			sq += exp((2.0*coordinates[mol].en_and_pr.energy - coordinates[mol].ex_field_coeff.energy)*beta - s);
			mob_histogram[mol] = sq;
		}

		for(int mol = 0; mol < nPart; mol++)
		{
			mob_histogram[mol] /= sq;
		}

		Rp = mob_histogram[nPart-1]*RanGen.Random();
		if(Rp < mob_histogram[0]) {trialPart = 0;}
			else
	        {
						int range = nPart/2;
						int point = range;
						do {
							if(range > 1) {range = range/2;}
							if (Rp > mob_histogram[point])
							{
								point = point + range;
							}
							else
							{
								point = point - range;
							}
						} while(!(Rp >= mob_histogram[point-1] && Rp < mob_histogram[point]));
						trialPart = point;
					}
	}
	state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle
	results delta_EP;
	results old_EP;
	results new_EP;
	results old_EP_Part;
	results new_EP_Part;
	// Move or Rotate a molecule
	new_coordinates.x = Lx * RanGen.Random();
	new_coordinates.y = Ly * RanGen.Random();
	new_coordinates.phi = 360.0*RanGen.Random();
	new_coordinates.sin_phi = sin(new_coordinates.phi/180.0*PI);
	new_coordinates.cos_phi = cos(new_coordinates.phi/180.0*PI);
  new_coordinates.damping_coeff = damping_field(new_coordinates.x, Lx); // Lambda^1/2
  new_coordinates.ex_field_coeff = external_field(new_coordinates.x, Lx); // u_ext
  new_coordinates.stat_weight = weights_for_central_cell (new_coordinates.x, Lx);
	charges_coordinates(new_coordinates);

for (int l = 0; l < nPart; l++)
  {
    check_HC(coordinates[l], new_coordinates, Lx, Ly);
  }
  if (HC_radius)
  {
    dt = 0.0;
  }
  else
  {
    for (int l = 0; l < nPart; l++)
      {
        if (l == trialPart){continue;}
        old_EP = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly, beta, true);
        new_EP = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, true);
        delta_EP = (new_EP - old_EP)/2.0;
        coordinates[l].en_and_pr = coordinates[l].en_and_pr + delta_EP;
        new_coordinates.en_and_pr = new_coordinates.en_and_pr + delta_EP;
      }
		new_coordinates.en_and_pr = new_coordinates.en_and_pr + (new_coordinates.ex_field_coeff - coordinates[trialPart].ex_field_coeff);
  	coordinates[trialPart] = new_coordinates;

		double s = 0;
		for(int mol = 0; mol < nPart; mol++)
		{
//			s += 2.0*coordinates[mol].en_and_pr.energy*beta;
			s += (2.0*coordinates[mol].en_and_pr.energy - coordinates[mol].ex_field_coeff.energy)*beta;
		}
		s = s/nPart;

		double sq = 0;
		for(int mol = 0; mol < nPart; mol++)
		{
			sq += exp((2.0*coordinates[mol].en_and_pr.energy - coordinates[mol].ex_field_coeff.energy)*beta - s);
			mob_histogram[mol] = sq;
		}
		for(int mol = 0; mol < nPart; mol++)
		{
			mob_histogram[mol] /= sq;
		}

		Rp = mob_histogram[nPart-1]*RanGen.Random();
		if(Rp < mob_histogram[0]) {trialPart = 0;}
			else
	        {
						int range = nPart/2;
						int point = range;
						do {
							if(range > 1) {range = range/2;}
							if (Rp > mob_histogram[point])
							{
								point = point + range;
							}
							else
							{
								point = point - range;
							}
						} while(!(Rp >= mob_histogram[point-1] && Rp < mob_histogram[point]));
						trialPart = point;
					}
		dt = exp(-s)/sq;
  }
  HC_radius = false;
  return dt;
}
