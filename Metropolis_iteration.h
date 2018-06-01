void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{
    double old_energy = EN_AND_PR_counter.energy;

    int trialPart = RanGen.IRandom(0,(nPart-1));
    double angle;
    bool angle_change = false;

    state old_coordinates = coordinates[trialPart];
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

      for (int l = 0; l < nPart; l++)
      {
           if (l == trialPart){continue;}
           //cout << endl << endl << " FIRST" << endl;
           delta_EP_old = energies_and_forces(coordinates[trialPart], coordinates[l], Lx, Ly, beta);
           double dopoln = delta_EP_old.energy;
           //cout << "SECOND" << endl;
           delta_EP_old = energies_and_forces(coordinates[l], coordinates[trialPart], Lx, Ly, beta);
           //if (abs(dopoln - delta_EP_old.energy) > abs(0.1*delta_EP_old.energy)) {cout << "AHTUNG!!!" << " dopoln=" << dopoln << " delta=" << delta_EP_old.energy << endl << endl;}
           old_EP = old_EP + delta_EP_old;
           delta_EP_new = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta);
           new_EP = new_EP + delta_EP_new;
           //add_E[l] = delta_EP_new.energy - delta_EP_old.energy;
           //add_E[trialPart] += add_E[l];
      }

      coordinates[trialPart] = new_coordinates;
      PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

      double new_energy = EN_AND_PR_counter.energy;
      delta_EP.energy = new_energy - old_energy;

      if (abs((new_EP.energy-old_EP.energy) - delta_EP.energy) > abs(delta_EP.energy*0.01))
      {
        cout << "diff=" << (new_EP.energy-old_EP.energy) << " exact=" << delta_EP.energy << endl;
        cout << "angle-" << angle_change << " old_mol=" << old_coordinates.x << " and " << old_coordinates.y << " and " << old_coordinates.phi << " new_mol=" << new_coordinates.x << " and " << new_coordinates.y << " and " << new_coordinates.phi << endl;
      }

      if(RanGen.Random() < exp(-delta_EP.energy*beta))
      {
          coordinates[trialPart] = new_coordinates;     // Update position
          //for(int l = 0; l < nPart; l++){coordinates[l].energy = coordinates[l].energy + add_E[l];}
          //EN_AND_PR_counter = EN_AND_PR_counter + delta_EP;
          if (angle_change) {ACCEPTANCE_RATIO_r[1] += 1.0;} else {ACCEPTANCE_RATIO_m[1] += 1.0;}
      }
      else
      {
          coordinates[trialPart] = old_coordinates;
          EN_AND_PR_counter.energy = old_energy;
          if(angle_change)
          {ACCEPTANCE_RATIO_r[0] += 1.0;}
          else
          {ACCEPTANCE_RATIO_m[0] += 1.0;}
      }
}
