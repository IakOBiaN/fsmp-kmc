#include <cassert>

void read_forcefield (const char * filename, vector <vector <vector <double> > > &forcefield, double &min_dist, double &max_dist, double &dr, double &da)
{
    //Create an input file stream
		FILE *pf = fopen (filename,"r");
		assert (pf != 0);
		assert (ferror(pf) == 0);

    vector<double> distance;
    vector<double> angle_1;
    vector<double> angle_2;
    vector<double> energy;

	  double line[4];
		while( 4 == fscanf (pf, "%lf %lf %lf %lf", line, line+1, line+2, line+3 ))
    {
				assert (ferror(pf) == 0);
				assert (feof(pf) == 0);
        distance.push_back(line[0]);
        angle_1.push_back(line[1]);
        angle_2.push_back(line[2]);
        energy.push_back(line[3]);
    }
		assert (feof(pf) != 0);

    //Close the file stream
		fclose(pf);

    double steps = distance.size();
    double angle_step_2;

    min_dist = distance[0];
		max_dist = distance[steps-1];

    for (int i=1;i<steps;i++)
    {
        if ((distance[i]-distance[0]) > 1e-5)
        {
           dr = distance[i]-distance[0];
           break;
        }
    }
    for (int i=1;i<steps;i++)
    {
        if ((angle_1[i]-angle_1[0]) > 1e-5)
        {
           da = angle_1[i]-angle_1[0];
           break;
        }
    }
    for (int i=1;i<steps;i++)
    {
        if ((angle_2[i]-angle_2[0]) > 1e-5)
        {
           angle_step_2 = angle_2[i]-angle_2[0];
           break;
        }
    }

    long int dist,ang1,ang2;
    for (long int i = 0; i < steps; i++) {
        dist = (int)(((distance[i]-distance[0])/dr)+0.5);
        ang1 = (int)(((angle_1[i]-angle_1[0])/da)+0.5);
        ang2 = (int)(((angle_2[i]-angle_2[0])/angle_step_2)+0.5);
        // 1 kcal = 4184 J/mol
       // forcefield[dist][ang1][ang2] = energy[i]*4184.0;
       forcefield[dist][ang1][ang2] = energy[i];
    }

			int a1, a2;
			for (int i = 0; i <= dist; i++)
			{
				for (int j = 0; j <= ang1; j++)
				{
					a1 = j + 180;
					if (a1 >= 360) {a1 -= 360;}
					for(int k = 0; k <= ang2; k++)
					{
						a2 = k + 180;
						if (a2 >= 360) {a2 -= 360;}
						if (abs(forcefield[i][j][k]-forcefield[i][a2][a1]) > 0.05*abs(forcefield[i][a2][a1]))
						{
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
