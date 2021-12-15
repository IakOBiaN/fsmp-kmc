using namespace std;

double density_to_Ly(double nPart, double state_dens)
{
	return sqrt((1.0e+26)*nPart/(state_dens*N_a*2.0/sqrt(3)));
}

double density_to_Ly_SF(double nPart, double state_dens)
{
	return sqrt((1.0e+26)*nPart/(state_dens*N_a*sqrt(3)));
}
