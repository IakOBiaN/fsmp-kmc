void weighted_averages(vector <state> &coordinates, int &nPart, double &Lx, double &Ly)
{
	nPart_in_central_cell = 0;
	for (int mol = 0; mol < nPart; mol++)
	{
		EN_AND_PR_counter.energy += coordinates[mol].en_and_pr.energy*coordinates[mol].stat_weight;
		EN_AND_PR_counter.p_X += coordinates[mol].en_and_pr.p_X*coordinates[mol].stat_weight;
		EN_AND_PR_counter.p_Y += coordinates[mol].en_and_pr.p_Y*coordinates[mol].stat_weight;
		nPart_in_central_cell += coordinates[mol].stat_weight;
	}
	density = (1.0e+26)*nPart_in_central_cell/(Lx/4.0*Ly)/N_a;
}
