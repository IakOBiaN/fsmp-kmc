void Metropolis_iteration(int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates)
{

    //vector <double> add_E(nPart);
    //for (int l = 0; l < nPart; l++){add_E[l] = 0;}

    int trialPart = RanGen.IRandom(0,(nPart-1));
    bool angle_change = false;

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    results delta_EP, old_EP, delta_EP_old, new_EP, delta_EP_new;

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
       new_coordinates.phi = coordinates[trialPart].phi + delta_angle*(2.0 * RanGen.Random() - 1.0);
       angle_change = true;
      }

    for (int l = 0; l < nPart; l++)
      {
           if (l == trialPart){continue;}
           delta_EP_old = energies_and_forces(coordinates[l], coordinates[trialPart], Lx, Ly, beta);
           old_EP = old_EP + delta_EP_old;
           delta_EP_new = energies_and_forces(coordinates[l], new_coordinates, Lx, Ly, beta);
           new_EP = new_EP + delta_EP_new;
           //add_E[l] = delta_EP_new.energy - delta_EP_old.energy;
           //add_E[trialPart] += add_E[l];

      }
      delta_EP = new_EP - old_EP;

      if(RanGen.Random() < exp(-delta_EP.energy))
      {
       // Accept rotation move
       coordinates[trialPart] = new_coordinates;     // Update angles
       //for(int l = 0; l < nPart; l++){coordinates[l].energy = coordinates[l].energy + add_E[l];}
       EN_AND_PR_counter = EN_AND_PR_counter + delta_EP;
       if (angle_change) {ACCEPTANCE_RATIO[1] += 1.0;}
      }
      else
        {if(angle_change) {ACCEPTANCE_RATIO[0] +=1.0;}}
}
