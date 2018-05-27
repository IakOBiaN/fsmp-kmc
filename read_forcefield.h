#include <cassert>

void read_forcefield (vector <vector <vector <double> > > &TMA_forcefield, double &min_dist, double &max_dist, double &dr, double &da)
{
    //Create an input file stream
		FILE *pf = fopen ("forcefield.dat","r");
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
    for (long int i = 0; i < steps; i++){
        dist = (distance[i]-distance[0])/dr;
        ang1 = (angle_1[i]-angle_1[0])/da;
        ang2 = (angle_2[i]-angle_2[0])/angle_step_2;
        // 1 kcal = 4184 J/mol
        TMA_forcefield.at(dist).at(ang1).at(ang2) = energy[i]*4184.0;
    }
}
