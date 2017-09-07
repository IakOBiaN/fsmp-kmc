void replace_the_trialParticle_and_update_energies(int &nPart, int &trialPart, double &Rc, double &Rc2, double &Lx, double &Ly,
                                                  double &beta, const double &A, double &C_q, vector <state> &coordinates)
{
 state new_coordinates; // Make a clone of trail particle

 // A random position is chosen uniformly over the whole volume of the system
 // New coordinates of the trial molecule
 // PBC is automatically satisfied
 new_coordinates.x = Lx * RanGen.Random();
 new_coordinates.y = Ly * RanGen.Random();
 new_coordinates.tetta = coordinates[trialPart].tetta + (3.141592653589/36)*(2.0 * RanGen.Random() - 1.0);
 new_coordinates.phi = coordinates[trialPart].phi + (3.141592653589/36)*(2.0 * RanGen.Random() - 1.0);
 new_coordinates.energy = 0;

 coordinates[trialPart].energy = 0;

 for (int l = 0; l < nPart; l++)
 {
    double g;
    if (l == trialPart){continue;}
    g = Inter_potential(coordinates[l], coordinates[trialPart], Rc, Rc2, Lx, Ly, A, C_q, beta);
    coordinates[l].energy -= g;
    g = Inter_potential(coordinates[l], new_coordinates, Rc, Rc2, Lx, Ly, A, C_q, beta);
    coordinates[l].energy += g;
    coordinates[trialPart].energy += g;
 }

 coordinates[trialPart].x = new_coordinates.x;
 coordinates[trialPart].y = new_coordinates.y;
 coordinates[trialPart].tetta = new_coordinates.tetta;
 coordinates[trialPart].phi = new_coordinates.phi;
 coordinates[trialPart].energy += (-334.4/((1/beta)*36.4))*abs(sin(coordinates[trialPart].tetta));
}
