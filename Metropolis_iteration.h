void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{
//    double CC_max_coord = (Lx/16.0)*9.0, CC_min_coord = (Lx/16.0)*7.0;

    int trialPart = RanGen.IRandom(0,(nPart-1));
    double angle;
    bool angle_change = false;

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    results delta_EP;
    results old_EP;
    results new_EP;
		results old_EP_Part;
		results new_EP_Part;

		if(RanGen.Random() < 0.5) // Move or Rotate a molecule
      {
       new_coordinates.x = coordinates[trialPart].x + (2 * delta * RanGen.Random() - delta); // random(-delta; delta)
       new_coordinates.y = coordinates[trialPart].y + (2 * delta * RanGen.Random() - delta);
       // Apply periodic boundary conditions
       new_coordinates.x = PBC2D(Lx, new_coordinates.x);
       new_coordinates.y = PBC2D(Ly, new_coordinates.y);
			 new_coordinates.damping_coeff = damping_field(new_coordinates.x, Lx); // Lambda^1/2
			 new_coordinates.ex_field_coeff = external_field(new_coordinates.x, Lx); // u_ext
			 new_coordinates.stat_weight = weights_for_central_cell (new_coordinates.x, Lx);
      }
    else
      {
       angle = coordinates[trialPart].phi + delta_angle*(2.0 * RanGen.Random() - 1.0);
       if (angle < 0) {angle += 360.0;}
       if (angle > 360) {angle -= 360.0;}
       new_coordinates.phi = angle;
       new_coordinates.sin_phi = sin(new_coordinates.phi/180.0*PI);
       new_coordinates.cos_phi = cos(new_coordinates.phi/180.0*PI);
       angle_change = true;
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
            //Choose exact or numerical energy and pressure calculation
            old_EP_Part = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta, true);
            new_EP_Part = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta, true);
            delta_EP = (new_EP_Part - old_EP_Part)/2.0;
            coordinates[l].en_and_pr = coordinates[l].en_and_pr + delta_EP;
    				new_coordinates.en_and_pr = new_coordinates.en_and_pr + delta_EP;
          }
					new_coordinates.en_and_pr = new_coordinates.en_and_pr + (new_coordinates.ex_field_coeff - coordinates[trialPart].ex_field_coeff);
          coordinates[trialPart] = new_coordinates;     // Update position
          if (angle_change) {ACCEPTANCE_RATIO_r[1] += 1.0;} else {ACCEPTANCE_RATIO_m[1] += 1.0;}
      }
			else
      {
				if(angle_change){ACCEPTANCE_RATIO_r[0] += 1.0;}
				else{ACCEPTANCE_RATIO_m[0] += 1.0;}
      }
      HC_radius = false;
}
