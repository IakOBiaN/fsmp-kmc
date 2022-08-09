void pressure_balance_area(double Energy, double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double S = Lx*Ly;
    dL = 0.001*min_dist;

		//pressure_balance for central cell
		if(press_X < press_Y)
		{
		 Ly_new = Ly + dL;
		 Lx_new = S/Ly_new;  // Area of the adlayer is kept constant
		}
		else
		{
		 Ly_new = Ly - dL;
		 Lx_new = S/Ly_new;
		}

		for(int i = 0; i < nPart; i++)
			 {
				 // Coordinates of the molecules change proportionally to
				 // the change of the box size
				 coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
				 coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
			 }
		Lx = Lx_new;
		Ly = Ly_new;

/*
			if (press_X < press_Y)
				{
					u_m += 0.15/beta;
				}
			else
				{
					u_m -= 0.15/beta;
				}
*/
       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

			 cout << "u: " << Energy/1000.0/nPart_in_central_cell << endl;
			 cout << "pX: " << (1.0e+23*nPart_in_central_cell/(Lx/4.0*Ly)/N_a/beta) + press_X*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " pY: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + press_Y*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " avg_p: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + (press_X + press_Y)*(1.0/(Lx/4.0)/Ly*1e23)/N_a/2.0 << endl;
			 cout << "Lx: " << Lx << " Ly: " << Ly << " density: " << nPart_in_central_cell*(1.0e+26)/(Lx/4.0*Ly)/N_a << " Um: " << u_m << endl;
}

void pressure_balance_ratio(double Energy, double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double Ly_Lx_ratio = Ly/Lx;
    dL = 0.001*min_dist;

		//pressure_balance for central cell
		if(press_X < press_Y)
		{
		 Ly_new = Ly + dL;
		 Lx_new = Ly_new/Ly_Lx_ratio;  // Area of the adlayer is kept constant
		}
		else
		{
		 Ly_new = Ly - dL;
		 Lx_new = Ly_new/Ly_Lx_ratio;
		}

		for(int i = 0; i < nPart; i++)
			 {
				 // Coordinates of the molecules change proportionally to
				 // the change of the box size
				 coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
				 coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
			 }
		Lx = Lx_new;
		Ly = Ly_new;

/*
			if (press_X < press_Y)
				{
					u_m += 0.15/beta;
				}
			else
				{
					u_m -= 0.15/beta;
				}
*/
       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

			 cout << "u: " << Energy/1000.0/nPart_in_central_cell << endl;
			 cout << "pX: " << (1.0e+23*nPart_in_central_cell/(Lx/4.0*Ly)/N_a/beta) + press_X*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " pY: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + press_Y*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " avg_p: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + (press_X + press_Y)*(1.0/(Lx/4.0)/Ly*1e23)/N_a/2.0 << endl;
			 cout << "Lx: " << Lx << " Ly: " << Ly << " density: " << nPart_in_central_cell*(1.0e+26)/(Lx/4.0*Ly)/N_a << " Um: " << u_m << endl;
}
