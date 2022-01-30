void zero_pressure_balance(double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
	double ratio_x_to_y = 2.0/sqrt(3); // for HON
	double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
	dL = 0.005;

	if(nPart/beta+(press_X + press_Y)/2.0 < 0)
		{
			Lx_new = Lx - dL;
		}
	else
		{
			Lx_new = Lx + dL;
		}
	Ly_new = Lx_new/ratio_x_to_y;

	// Coordinates of the molecules change proportionally to
	// the change of the box size
	for(int i = 0; i < nPart; i++)
		{
			coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
			coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
		}

	// Change the size of the simulation cell
	Lx = Lx_new;
	Ly = Ly_new;
	density = (1.0e+26)*nPart/(Lx*Ly)/N_a;

	// Recalculate energies after compressing or expanding the box
	PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

	cout << "pX: " << (1.0e+23*nPart/(Lx*Ly)/N_a/beta) + press_X*(1.0/Lx/Ly*1e23)/N_a << " pY: " << (1.0e+23*nPart/(Lx*Ly)/N_a/beta) + press_Y*(1.0/Lx/Ly*1e23)/N_a << " avg_p: " << (1.0e+23*nPart/(Lx*Ly)/N_a/beta) + (press_X + press_Y)*(1.0/Lx/Ly*1e23)/N_a/2.0 << endl;
	cout << "Lx: " << Lx << " Ly: " << Ly << " density: " << density << endl;
}
