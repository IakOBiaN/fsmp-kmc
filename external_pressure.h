using namespace std;

double external_pressure(double x, double Lx, double potential)
{
 double p = 0;
 double v = 0;
 double v0 = 2.0;
 double c = potential;
 if(x <= 0.5*Lx){v = 8.0*(0.5*Lx - x)/Lx;}
	else if(x > 0.5*Lx){v = 8.0*(x - 0.5*Lx)/Lx;}
		if(v > v0){p = 4*c*(v - 2.0)/Lx;}
			else if(v <= v0) {p = 0;}
 return 0;
}
