void weighted_averages_in_central_cell(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	nPart_in_central_cell = 0;
	for (int mol = 0; mol < nPart; mol++)
	{
		EN_AND_PR_counter.energy += coordinates[mol].en_and_pr.energy*coordinates[mol].stat_weight;
		EN_AND_PR_counter.p_X += coordinates[mol].en_and_pr.p_X*coordinates[mol].stat_weight;
		EN_AND_PR_counter.p_Y += coordinates[mol].en_and_pr.p_Y*coordinates[mol].stat_weight;
		nPart_in_central_cell += coordinates[mol].stat_weight;
	}
}

void weighted_averages_in_gas(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	EN_AND_PR_counter.energy = 0;
	EN_AND_PR_counter.p_X = 0;
	EN_AND_PR_counter.p_Y = 0;
	nPart_in_gas = 0;
	for (int mol = 0; mol < nPart; mol++)
	{
		double abs_x;
		if (coordinates[mol].x > Lx/2.0){abs_x = coordinates[mol].x - Lx/2.0;} else {abs_x = Lx/2.0 - coordinates[mol].x;}
		double ksi = 32.0*abs_x/Lx - 8.0;
		if(ksi > 4.0)
			{
				EN_AND_PR_counter.energy += coordinates[mol].en_and_pr.energy;
				EN_AND_PR_counter.p_X += coordinates[mol].en_and_pr.p_X;
				EN_AND_PR_counter.p_Y += coordinates[mol].en_and_pr.p_Y;
				nPart_in_gas += 1.0;
			}
	}
}
