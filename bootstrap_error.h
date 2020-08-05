using namespace std;

double bootstrap_error_calculation(vector <double> data, double num_elements)
{
	int n = 500;
	double error_estimation;
	vector <double> avrg(n);
	for(int i = 0; i < n; i++)
		{
			double block_avrg = 0;
			for(int j = 1; j <= num_elements; j++)
				{
					int element = RanGen.IRandom(1,num_elements);
					block_avrg += data[element];
				}
			avrg[i] = block_avrg/num_elements;
		}

	double avrg_value = 0;
	double avrg_sqr_value = 0;
	for(int k = 0; k < n; k++)
		{
			avrg_value += avrg[k];
			avrg_sqr_value += avrg[k]*avrg[k];
		}
	avrg_value /= n;
	avrg_sqr_value /= n;
	error_estimation = sqrt(avrg_sqr_value - avrg_value*avrg_value);
	return error_estimation;
}
