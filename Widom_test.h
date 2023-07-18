void Widom_test(int &nPart, vector <state> &coordinates, double &Lx, double &Ly, double &beta, int &N_test, double &e_test)
{
  state test_particle;
  double test_particle_energy = 0;

	if (RanGen.Random() < 0.5){test_particle.x = (7.0 + RanGen.Random())*Lx/8.0;}
	else {test_particle.x = RanGen.Random()*Lx/8.0;}
  //test_particle.x = Lx * RanGen.Random();
  test_particle.y = Ly * RanGen.Random();
  test_particle.phi = 360.0 * RanGen.Random();
  test_particle.sin_phi = sin(test_particle.phi/180.0*PI);
  test_particle.cos_phi = cos(test_particle.phi/180.0*PI);
	test_particle.damping_coeff = damping_field(test_particle.x, Lx); // Lambda^1/2
	test_particle.ex_field_coeff = external_field(test_particle.x, Lx); // u_ext
	test_particle.stat_weight = weights_for_central_cell(test_particle.x, Lx);
	for (int l = 0; l < nPart; l++){test_particle_energy = test_particle_energy +  energies_and_forces(test_particle, coordinates[l], Lx, Ly, beta, false).energy;}
	test_particle_energy = test_particle_energy + test_particle.ex_field_coeff.energy;
	e_test += exp(-test_particle_energy*beta);
  N_test++;
  HC_radius = false;
}
