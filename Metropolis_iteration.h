void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{
    int trialPart = RanGen.IRandom(0,(nPart-1));
    double angle;
    bool angle_change = false;

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    results delta_EP;
    results old_EP;
    results delta_EP_old;
    results new_EP;
    results delta_EP_new;
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

      //#pragma omp parallel for reduction(+:old_EP) reduction(+:new_EP)
      for (int l = 0; l < nPart; l++)
      {
           if (l == trialPart){continue;}
           old_EP = old_EP +  energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly,beta);
           //old_EP = old_EP + delta_EP_old;
           new_EP = new_EP + energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta);
           //new_EP = new_EP + delta_EP_new;
      }

      delta_EP = new_EP - old_EP;
      if(RanGen.Random() < exp(-delta_EP.energy))
      {
          coordinates[trialPart] = new_coordinates;     // Update position
          EN_AND_PR_counter = EN_AND_PR_counter + delta_EP;
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
