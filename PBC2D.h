using namespace std;

double PBC2D (double &L, double d)
{
 if(d > L){d = d - int(d/L)*L;}
  else if(d < 0){d = d + L;}
 return d;
}
