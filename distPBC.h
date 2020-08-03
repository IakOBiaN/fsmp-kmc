using namespace std;

// L is the box size
// d is a particle-particle distance along the X, Y axis
double distPBC (double &L, double d)
{
 d = fmod(d, L);
 if(d > (L - d)){d = L - d;}
 return d;
}
