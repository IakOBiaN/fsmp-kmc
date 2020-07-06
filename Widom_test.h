void Widom_test(int &nPart, vector <state> &coordinates, double &Lx, double &Ly, double &beta, int &N_test, double &e_test)
{
  state test_particle;
  double test_particle_energy = 0;

  test_particle.x = Lx * RanGen.Random();
  test_particle.y = Ly * RanGen.Random();
  test_particle.phi = 360.0 * RanGen.Random();
  test_particle.sin_phi = sin(test_particle.phi/180.0*PI);
  test_particle.cos_phi = cos(test_particle.phi/180.0*PI);
  charges_coordinates(test_particle);
  for (int l = 0; l < nPart; l++){test_particle_energy = test_particle_energy +  energies_and_forces_exact(test_particle, coordinates[l], Lx, Ly, beta).energy;}
//	for (int l = 0; l < nPart; l++){test_particle_energy = test_particle_energy +  energies_and_forces(test_particle, coordinates[l], Lx, Ly, beta, true).energy;}
  e_test += exp(-test_particle_energy*beta);
  N_test++;
}
