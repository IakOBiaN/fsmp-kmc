void Metropolis_iteration(int &nPart, double &Rc, double &Rc2, double &Lx, double &Ly,
                          double &beta, const double &A, double &C_q, vector <state> &coordinates)
{
    double delta = 0.5;
    double delta_angle = 3.141592653589/36; // Maximal rotation in rad

    vector <double> add_E(5000);
    for (int l = 0; l < nPart; l++){add_E[l] = 0;}

    int trialPart = RanGen.IRandom(0,(nPart-1));

    state new_coordinates = coordinates[trialPart]; // Make a clone of trail particle

    double deltaE = 0;
    double oldE = 0, delta_E_old = 0;
    double newE = 0, delta_E_new = 0;

    if(RanGen.Random() < 0.5) // Move or Rotate a molecule
      {
       new_coordinates.x = coordinates[trialPart].x + (2 * delta * RanGen.Random() - delta); // random(-delta; delta)
       new_coordinates.y = coordinates[trialPart].y + (2 * delta * RanGen.Random() - delta);

       // Apply periodic boundary conditions
       new_coordinates.x = PBC2D(Lx, new_coordinates.x);
       new_coordinates.y = PBC2D(Ly, new_coordinates.y);

       for (int l = 0; l < nPart; l++)
          {
            if (l == trialPart){continue;}
            delta_E_old = Inter_potential(coordinates[l], coordinates[trialPart], Rc, Rc2, Lx, Ly, A, C_q, beta);
            oldE += delta_E_old;
            delta_E_new = Inter_potential(coordinates[l], new_coordinates, Rc, Rc2, Lx, Ly, A, C_q, beta);
            newE += delta_E_new;
            add_E[l] = delta_E_new - delta_E_old;
            add_E[trialPart] += add_E[l];
          }
       deltaE = newE - oldE;

       if(RanGen.Random() < exp(-beta*deltaE))
         {
          // Accept displacement move
          coordinates[trialPart] = new_coordinates;     // Update coordinates
          for(int l = 0; l < nPart; l++){coordinates[l].energy = coordinates[l].energy + add_E[l];}
         }
      }
    else
      {
       new_coordinates.phi = coordinates[trialPart].phi + delta_angle*(2.0 * RanGen.Random() - 1.0);

       for (int l = 0; l < nPart; l++)
          {
            if (l == trialPart){continue;}
            delta_E_old = Inter_potential(coordinates[l], coordinates[trialPart], Rc, Rc2, Lx, Ly, A, C_q, beta);
            oldE += delta_E_old;
            delta_E_new = Inter_potential(coordinates[l], new_coordinates, Rc, Rc2, Lx, Ly, A, C_q, beta);
            newE += delta_E_new;
            add_E[l] = delta_E_new - delta_E_old;
            add_E[trialPart] += add_E[l];
          }
       deltaE = newE - oldE;

       if(RanGen.Random() < exp(-deltaE))
         {
          // Accept rotation move
          coordinates[trialPart] = new_coordinates;     // Update angles
          for(int l = 0; l < nPart; l++){coordinates[l].energy = coordinates[l].energy + add_E[l];}
         }
      }
}
