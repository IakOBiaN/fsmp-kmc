using namespace std;

double PBC2D (double &L, double d)
{
 if(d > L){d = d - L;}
  else if(d < 0){d = d + L;}
 return d;
}
