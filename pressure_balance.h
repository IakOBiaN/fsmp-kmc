void pressure_balance_area(double Energy, double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double S = Lx*Ly;
    dL = 0.0001*Ly;

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
    dL = 0.00005*Ly;

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
				 coordinates[i].damping_coeff = damping_field(coordinates[i].x, Lx_new); // Lambda^1/2
				 coordinates[i].ex_field_coeff = external_field(coordinates[i].x, Lx_new); // u_ext
				 coordinates[i].stat_weight = weights_for_central_cell (coordinates[i].x, Lx_new);
				 charges_coordinates (coordinates[i]);
			 }
		Lx = Lx_new;
		Ly = Ly_new;

       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

			 cout << "u: " << Energy/1000.0/nPart_in_central_cell << endl;
			 cout << "pX: " << (1.0e+23*nPart_in_central_cell/(Lx/4.0*Ly)/N_a/beta) + press_X*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " pY: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + press_Y*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " avg_p: " << (1.0e+23*nPart_in_central_cell/((Lx/4.0)*Ly)/N_a/beta) + (press_X + press_Y)*(1.0/(Lx/4.0)/Ly*1e23)/N_a/2.0 << endl;
			 cout << "Lx: " << Lx << " Ly: " << Ly << " density: " << nPart_in_central_cell*(1.0e+26)/(Lx/4.0*Ly)/N_a << " Um: " << u_m << endl;
}

void pressure_balance_ratio_analytical(double Energy, double press_X, double delta_p_over_interface, double gas_density, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double Ly_Lx_ratio = Ly/Lx;
    dL = 0.0001*Ly;

		//pressure_balance for central cell
		if(press_X*4.0/Lx < delta_p_over_interface)
		{
		 Lx_new = Lx - dL;
		 Ly_new = Lx_new*Ly_Lx_ratio;
		}
		else
		{
		 Lx_new = Lx + dL;
		 Ly_new = Lx_new*Ly_Lx_ratio;
		}

		for(int i = 0; i < nPart; i++)
			 {
				 // Coordinates of the molecules change proportionally to
				 // the change of the box size
				 coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
				 coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
				 coordinates[i].damping_coeff = damping_field(coordinates[i].x, Lx_new); // Lambda^1/2
				 coordinates[i].ex_field_coeff = external_field(coordinates[i].x, Lx_new); // u_ext
				 coordinates[i].stat_weight = weights_for_central_cell (coordinates[i].x, Lx_new);
				 charges_coordinates (coordinates[i]);
			 }
		Lx = Lx_new;
		Ly = Ly_new;

		// Recalculate energies after compressing or expanding the box
		PotentialEnergy(nPart, Lx, Ly, coordinates, beta);

			 cout << "u: " << Energy/1000.0/nPart_in_central_cell << endl;
			 cout << "pX: " << (1.0e+23*nPart_in_central_cell/(Lx/4.0*Ly)/N_a/beta) + press_X*(1.0/(Lx/4.0)/Ly*1e23)/N_a << " dp: " << R*temperature*gas_density/1000.0 + delta_p_over_interface*1e23/Ly/N_a << endl;
			 cout << "Lx: " << Lx << " Ly: " << Ly << " density: " << nPart_in_central_cell*(1.0e+26)/(Lx/4.0*Ly)/N_a << " Um: " << u_m << endl;
}

void um_tunning_to_zero_pressure(double Energy, double &u_m, double delta_p_over_interface, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
		double d_um = 500;
		//pressure_balance for central cell
		if(R*temperature*gas_density/1000.0 + delta_p_over_interface*1e23/Ly/N_a > 0)
		{
			u_m -= d_um;
		}
		else
		{
			u_m += d_um;
		}
		for(int i = 0; i < nPart; i++){coordinates[i].ex_field_coeff = external_field(coordinates[i].x, Lx);}


       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
			 cout << "u: " << Energy/1000.0/nPart_in_central_cell << endl;
			 cout << "Um: " << u_m << endl;
}
