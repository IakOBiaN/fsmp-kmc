void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates, double potential)
{
//    double CC_max_coord = (Lx/16.0)*9.0, CC_min_coord = (Lx/16.0)*7.0;

    int trialPart = RanGen.IRandom(0,(nPart-1));
    double angle;
    bool angle_change = false;

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    results delta_EP;
		results delta_EP_cent;
    results old_EP;
    results new_EP;
		results old_EP_Part;
		results new_EP_Part;
		results old_EP_Part_cent;
		results new_EP_Part_cent;
    results old_EP_cent;
		results new_EP_cent;

		double delta_ext_press = 0;

    if(RanGen.Random() < 0.5) // Move or Rotate a molecule
      {
       new_coordinates.x = coordinates[trialPart].x + (2 * delta * RanGen.Random() - delta); // random(-delta; delta)
       new_coordinates.y = coordinates[trialPart].y + (2 * delta * RanGen.Random() - delta);
       // Apply periodic boundary conditions
       new_coordinates.x = PBC2D(Lx, new_coordinates.x);
       new_coordinates.y = PBC2D(Ly, new_coordinates.y);
			 //new_coordinates.cent = cent_potential(new_coordinates.x, Lx);
       charges_coordinates(new_coordinates);
      }
    else
      {
       angle = coordinates[trialPart].phi + delta_angle*(2.0 * RanGen.Random() - 1.0);
       if (angle < 0) {angle += 360.0;}
       if (angle > 360) {angle -= 360.0;}
       new_coordinates.phi = angle;
       new_coordinates.sin_phi = sin(new_coordinates.phi/180.0*PI);
       new_coordinates.cos_phi = cos(new_coordinates.phi/180.0*PI);
       charges_coordinates(new_coordinates);
       angle_change = true;
      }

    for (int l = 0; l < nPart; l++)
    {
         if (l == trialPart){continue;}
         old_EP_Part = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta,true);
         new_EP_Part = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta,true);
				 old_EP_Part_cent = old_EP_Part*coordinates[trialPart].cent * coordinates[l].cent;
				 new_EP_Part_cent = new_EP_Part*new_coordinates.cent*coordinates[l].cent;
				 old_EP = old_EP + old_EP_Part;
				 new_EP = new_EP + new_EP_Part;
				 //old_EP_cent = old_EP_cent + old_EP_Part_cent;
				 //new_EP_cent = new_EP_cent + new_EP_Part_cent;
		}
    delta_EP = new_EP - old_EP;
    //delta_EP_cent = new_EP_cent - old_EP_cent;
    //delta_EP.energy +=  (external_field(new_coordinates.x, Lx,potential) - external_field(coordinates[trialPart].x, Lx, potential));
    if(RanGen.Random() < exp(-delta_EP.energy*beta))
     {
	    EN_AND_PR_counter = EN_AND_PR_counter + delta_EP;
      EN_AND_PR_counter_central_cell = EN_AND_PR_counter;
			//EN_AND_PR_counter_central_cell = EN_AND_PR_counter_central_cell + delta_EP_cent;
			//ext_pressure += (external_pressure(new_coordinates.x, Lx,potential) - external_pressure(coordinates[trialPart].x, Lx, potential));
	     //cout << "Energy:" << EN_AND_PR_counter.energy << " P_x:" << EN_AND_PR_counter.p_X << " P_y:" << EN_AND_PR_counter.p_Y << endl;
	     //cout << "nEnergy:" << delta_EP.energy << " P_x:" << delta_EP.p_X << " P_y:" << delta_EP.p_Y << endl;
	    coordinates[trialPart] = new_coordinates;     // Update position
	    if (angle_change) {ACCEPTANCE_RATIO_r[1] += 1.0;} else {ACCEPTANCE_RATIO_m[1] += 1.0;}
     }
    else
     {
      if(angle_change){ACCEPTANCE_RATIO_r[0] += 1.0;}
       else {ACCEPTANCE_RATIO_m[0] += 1.0;}
     }
}
