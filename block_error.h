using namespace std;

double block_error_calculation(vector <double> data, double num_elements)
{
	int n = 10;
	int period = num_elements/n;
	double error_estimation;
	vector <double> avrg(n);
	vector <double> avrg_sqr(n);
	double block_avrg = 0;
	double block_avrg_sqr = 0;
	int j = 0;
	for(int i = 1; i <= num_elements; i++)
			{
				if (i%period == 0){j += 1;}
				avrg[j] += data[i];
				avrg_sqr[j] += data[i]*data[i];
			}
	for(int k = 0; k < n; k++)
			{
				block_avrg += avrg[k]/period;
				block_avrg_sqr += avrg_sqr[k]/period;
			}
	block_avrg /= n;
	block_avrg_sqr /= n;

	error_estimation = sqrt((block_avrg_sqr - block_avrg*block_avrg)/(n-1));
	return error_estimation;
}
