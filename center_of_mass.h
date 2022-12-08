double center_of_mass(int &nPart, vector <state> &coordinates)
{
	// Calculate x-coordinate of mass center of the system
	double xc = 0;
	for (int i = 0; i < nPart; i++)
		{
			xc += coordinates[i].x;
		}
	xc = xc/nPart;
	return xc;
}
