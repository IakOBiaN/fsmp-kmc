using namespace std;

double external_field(double x, double Lx)
{
 double f0 = 0;
 //double f1 = 0;
 double v = 0;
 double v0 = 2.0;
 //double v1 = 3.0;
 double c = 100000;
 if(x <= 0.5*Lx){v = 8.0*(0.5*Lx - x)/Lx;}
  else if(x > 0.5*Lx){v = 8.0*(x - 0.5*Lx)/Lx;}
        if(v > v0){f0 = c*pow((v - 2.0), 2)/4.0;}
         else if(v <= v0) {f0 = 0;}
        //if(v > v1){f1 = -v + 3.0 - log(4.0 - v);}
         //else if(v <= v1) {f1 = 0;}
 //return (f0 + f1);
 return f0;
}
