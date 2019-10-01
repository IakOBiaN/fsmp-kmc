void Widom_test(int &nPart, vector <state> &coordinates, double &Lx, double &Ly, double &beta, int &N_test, double &e_test)
{
  state test_particle;
  double test_particle_energy = 0;

  if((N_test%2) == 0)
    {
      test_particle.x = 0.25*Lx * RanGen.Random();
    }
    else {test_particle.x = Lx * (1.0 - 0.25 * RanGen.Random());}
  test_particle.y = Ly * delta * RanGen.Random();
  test_particle.phi = 360.0 * RanGen.Random();

  for (int l = 0; l < nPart; l++){test_particle_energy = test_particle_energy +  energies_and_forces(test_particle, coordinates[l], Lx, Ly,beta,false).energy;}

  test_particle_energy += external_field(test_particle.x, Lx);
  e_test += exp(-test_particle_energy*beta);
  N_test ++;
}
