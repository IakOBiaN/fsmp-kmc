#include <cassert>

void read_forcefield (const char * filename, vector <vector <vector <double> > > &forcefield, double &min_dist, double &max_dist, double &dr, double &da)
{
    //Create an input file stream
		FILE *pf = fopen (filename,"r");
		assert (pf != 0);
		assert (ferror(pf) == 0);

	  double line[4];
		double first_dist = -1;
		double first_a1;
		double first_a2;
		double now_dist = -100;
		double now_a1 = -100;
		double now_a2 = -100;
		int i = 0;
		int j = 0;
		int k = -1;
		int counter_dist = 0;
		int counter_da = 0;
		int counter_da2 = 0;
		while( 4 == fscanf (pf, "%lf %lf %lf %lf", line, line+1, line+2, line+3 ))
    {
				assert (ferror(pf) == 0);
				assert (feof(pf) == 0);
				if (first_dist == -1)
				{
					first_dist = line[0];
					first_a1 = line[1];
					first_a2 = line[2];

					now_a1 = line[1];
					now_dist = line[0];
				}

				k++;
				if (now_a1 != line[1])
				{
					da = line[1] - now_a1;
					k = 0;
					j++;
					if (now_dist != line[0])
					{
						cout << "Now r:" << line[0] << endl;
						dr = line[0] - now_dist;
						j = 0;
						i++;
					}
					now_dist = line[0];
				}
				now_a1 = line[1];

				forcefield[i][j][k] = line[3] * 4184.0;

    }
		assert (feof(pf) != 0);

    //Close the file stream
		fclose(pf);

		min_dist = first_dist;
		max_dist = now_dist;

		cout << "min:" << min_dist << " max:" << max_dist << endl;
		cout << "dr:" << dr << " da:" << da << endl;
		cout << "max_ijk:" << i << " " << j << " " << k << endl;

			int dist = i;
			int ang1 = j;
			int ang2 = k;
			int a1, a2;
			for (int i = 0; i <= dist; i++)
			{
				for (int j = 0; j <= ang1; j++)
				{
					a1 = j + int(180 / da + 0.5);
					if (a1 >= int(360 / da + 0.5)) {a1 -= int(360 / da + 0.5);}
					for(int k = 0; k <= ang2; k++)
					{
						a2 = k + int(180 / da + 0.5);
						if (a2 >= int(360 / da + 0.5)) {a2 -= int(360 / da + 0.5);}
						// Test of the potential for symmetry
						//if (abs(forcefield[i][j][k]-forcefield[i][a2][a1]) > 0.15*abs(forcefield[i][a2][a1]))
						if (abs(forcefield[i][j][k]-forcefield[i][a2][a1]) > 1000 && abs(forcefield[i][j][k]) < E_INF * R * 1000.0)
								{
									cout << "distance: " << min_dist + i * dr << " angle1: " << j * da << " angle2: " << k * da << endl;
									cout << "Potential error. Energy=" << forcefield[i][j][k] << " must be close energy=" << forcefield[i][a2][a1] << endl;
								}
						if (forcefield[i][j][k] < forcefield[i][a2][a1])
						{
							forcefield[i][j][k] = forcefield[i][a2][a1];
						}
						else
						{
							forcefield[i][a2][a1] = forcefield[i][j][k];
						}
					}
				}
			}
}
