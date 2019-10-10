void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{
    double CC_max_coord = (Lx/16.0)*9.0, CC_min_coord = (Lx/16.0)*7.0;

    int trialPart = RanGen.IRandom(0,(nPart-1));
    double angle;
    bool angle_change = false;

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    results delta_EP;
    results old_EP;
    results new_EP_central;
    results old_EP_central;
    results delta_EP_old;
    results delta_EP_new;
    results EP_with_central;
    bool out_to_in,in_to_out;
    results new_EP;
    if(RanGen.Random() < 0.5) // Move or Rotate a molecule
      {
       new_coordinates.x = coordinates[trialPart].x + (2 * delta * RanGen.Random() - delta); // random(-delta; delta)
       new_coordinates.y = coordinates[trialPart].y + (2 * delta * RanGen.Random() - delta);
       // Apply periodic boundary conditions
       new_coordinates.x = PBC2D(Lx, new_coordinates.x);
       new_coordinates.y = PBC2D(Ly, new_coordinates.y);
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

      out_to_in = false;
      if ((coordinates[trialPart].x < CC_min_coord || coordinates[trialPart].x > CC_max_coord) && (new_coordinates.x > CC_min_coord && new_coordinates.x < CC_max_coord))
      {
        out_to_in = true;
      }
      in_to_out = false;
      if ((coordinates[trialPart].x > CC_min_coord && coordinates[trialPart].x < CC_max_coord) && (new_coordinates.x < CC_min_coord || new_coordinates.x > CC_max_coord))
      {
        in_to_out = true;
      }
      //#pragma omp parallel for reduction(+:old_EP) reduction(+:new_EP)
      for (int l = 0; l < nPart; l++)
      {
           if (l == trialPart){continue;}

           delta_EP_old = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta,true);
           delta_EP_new = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta,true);
					 //Choose exact or numerical energy and pressure calculation
          if (out_to_in)
          {
            if (coordinates[l].x > CC_min_coord && coordinates[l].x < CC_max_coord)
            {
                EP_with_central = EP_with_central + delta_EP_old;
            }
          }
          if (in_to_out)
          {
            if (coordinates[l].x > CC_min_coord && coordinates[l].x < CC_max_coord)
            {
                EP_with_central = EP_with_central + delta_EP_new;
            }
          }
            if (coordinates[l].x > CC_min_coord && coordinates[l].x < CC_max_coord)
            {
                old_EP_central = old_EP_central + delta_EP_old;
                new_EP_central = new_EP_central + delta_EP_new;
            }


           old_EP = old_EP +  delta_EP_old;
           new_EP = new_EP + delta_EP_new;
      }

      delta_EP = new_EP - old_EP;
      delta_EP.energy +=  (external_field(new_coordinates.x, Lx) - external_field(coordinates[trialPart].x, Lx));
      //cout << "Metropolis_energy: " << delta_EP.energy*beta << endl;
      if(RanGen.Random() < exp(-delta_EP.energy*beta))
      {
          if ((coordinates[trialPart].x > CC_min_coord && coordinates[trialPart].x < CC_max_coord) && (new_coordinates.x > CC_min_coord && new_coordinates.x < CC_max_coord))
              {
                EN_AND_PR_counter = EN_AND_PR_counter + delta_EP;
                //cout << "Energy:" << EN_AND_PR_counter.energy << " P_x:" << EN_AND_PR_counter.p_X << " P_y:" << EN_AND_PR_counter.p_Y << endl;
                //cout << "nEnergy:" << delta_EP.energy << " P_x:" << delta_EP.p_X << " P_y:" << delta_EP.p_Y << endl;
              }
          else if (out_to_in)
              {
                EN_AND_PR_counter = EN_AND_PR_counter + new_EP - EP_with_central;
              }
          else if (in_to_out)
              {
                EN_AND_PR_counter = EN_AND_PR_counter + EP_with_central - old_EP;
              }
          else
              {
                EN_AND_PR_counter = EN_AND_PR_counter + new_EP_central - old_EP_central;
              }
          coordinates[trialPart] = new_coordinates;     // Update position
          if (angle_change) {ACCEPTANCE_RATIO_r[1] += 1.0;} else {ACCEPTANCE_RATIO_m[1] += 1.0;}
      }
      else
      {
          if(angle_change)
          {ACCEPTANCE_RATIO_r[0] += 1.0;}
          else
          {ACCEPTANCE_RATIO_m[0] += 1.0;}
      }
}
