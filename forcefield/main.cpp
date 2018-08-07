#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <float.h>

using namespace std;

class pressure {
public:
double X_LJ;
double X_QQ;
double Y_LJ;
double Y_QQ;
pressure operator+(pressure& b) {
         pressure press;
         press.X_LJ = this->X_LJ + b.X_LJ;
         press.X_QQ = this->X_QQ + b.X_QQ;
         press.Y_LJ = this->Y_LJ + b.Y_LJ;
         press.Y_QQ = this->Y_QQ + b.Y_QQ;
         return press;
      }
pressure operator-(pressure& b) {
         pressure press;
         press.X_LJ = this->X_LJ - b.X_LJ;
         press.X_QQ = this->X_QQ - b.X_QQ;
         press.Y_LJ = this->Y_LJ - b.Y_LJ;
         press.Y_QQ = this->Y_QQ - b.Y_QQ;
         return press;
      }
};

class results {
public:
double energy;
double energy_QQ;
pressure p;
results();      //constructor
results operator+(results& b) {
         results res;
         res.energy = this->energy + b.energy;
         res.energy_QQ = this->energy_QQ + b.energy_QQ;
         res.p = this->p + b.p;
         return res;
      }
results operator-(results& b) {
         results res;
         res.energy = this->energy - b.energy;
         res.energy_QQ = this->energy_QQ - b.energy_QQ;
         res.p = this->p - b.p;
         return res;
      }
};

//constructor
results::results(void) {
   energy = 0;
   energy_QQ = 0;
   p.X_LJ = 0;
   p.X_QQ = 0;
   p.Y_LJ = 0;
   p.Y_QQ = 0;
}

class state {
public:
double x;
double y;
double phi;
double sin_phi;
double cos_phi;
double energy;
double mob;
};
double PI = 3.1415926535;

double R = 8.3144598;
double eps = 0.502e-21;                              // LJ energy for nitrogen in J
double N_a = 6.02214e+23;
double sigma = 3.318e-10;                                // in meters
double k_B = 1.38e-23;
double Rc = 20.0;                              // Cut-off radius in sigma
double Rc2 = Rc*Rc;
double Qn2 = -4.453e-40;             // Quadrupole moment of N2 molecule
const double eps0 = 8.85418781762e-12;                     // The permittivity of free space in C2 m-2 N-1
const double A = 1.0/(4.0*3.1415926535*eps0);    // Coulomb's constant
double C_q = A*(3.0/4.0)*pow(Qn2,2);
double dn2 = 0.33092224232;               // Distance between nitrogen atoms in sigma units
double dq1 = 0.25527426160;               // Distance between "+" charge and center of quadrupole in sigma units
double dq2 = 0.31464737794;               // Distance between "-" charge and center of quadrupole in sigma units
const double qe = 1.6021766208e-19;                   // The charge of one electron in C
double q = 0.373*qe;                            // Charge of the quadrupole points in C

#include "energies_and_forces.h"

int main()
{
    state molA, molB;
    state trial_Part;
    molA.x = 0.0;
    molA.y = 0.0;
    molB.y = 0.0;
    results pair_energy;
    results trial_energy_4, trial_energy_3, trial_energy_2, trial_energy_1, trial_energy1, trial_energy2, trial_energy3, trial_energy4;
    double pair_force_LJ, pair_force_QQ;

    stringstream name;
    name <<  "forcefield.dat";
    ofstream fileOutput(name.str().c_str(), ios_base::trunc);

    stringstream name1;
    name1 <<  "force_LJ.dat";
    ofstream fileOutput1(name1.str().c_str(), ios_base::trunc);

    stringstream name2;
    name2 <<  "force_QQ.dat";
    ofstream fileOutput2(name2.str().c_str(), ios_base::trunc);

    double min_dist = 2.0;

    for (int i=0; i < 131; i++) // distance
    {
        cout << "i: " << i << endl;
        molB.x = (min_dist + i*0.1)*(1e-10); // in meters
        for (int j=0; j < 361; j+=1.0)
        {
            molA.phi = j;
            molA.sin_phi = sin(molA.phi/180.0*PI);
            molA.cos_phi = cos(molA.phi/180.0*PI);
            for (int k=0; k < 361; k+=1.0)
            {
                double dFi_dr;
                double delta_r = 0.1/(1e+10);


                molB.phi = k;
                molB.sin_phi = sin(molB.phi/180.0*PI);
                molB.cos_phi = cos(molB.phi/180.0*PI);
                pair_energy = energies_and_forces(molA, molB, 1000.0, 1000.0);

                double r_0 = sqrt(pow((molA.x - molB.x), 2) + pow((molA.y - molB.y), 2));
                fileOutput << r_0*(1e+10) << " " << molA.phi << " " << molB.phi << " " << pair_energy.energy + pair_energy.energy_QQ <<  endl;

                trial_Part = molB;
                trial_Part.x = molB.x - 4.0*delta_r;
                trial_energy_4 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x - 3.0*delta_r;
                trial_energy_3 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x - 2.0*delta_r;
                trial_energy_2 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x - 1.0*delta_r;
                trial_energy_1 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x + 1.0*delta_r;
                trial_energy1 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x + 2.0*delta_r;
                trial_energy2 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x + 3.0*delta_r;
                trial_energy3 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);
                trial_Part.x = molB.x + 4.0*delta_r;
                trial_energy4 = energies_and_forces(molA, trial_Part, 1000.0, 1000.0);

                dFi_dr = (1.0/280.0*trial_energy_4.energy - 4.0/105.0*trial_energy_3.energy + 1.0/5.0*trial_energy_2.energy - 4.0/5.0*trial_energy_1.energy + 4.0/5.0*trial_energy1.energy - 1.0/5.0*trial_energy2.energy + 4.0/105.0*trial_energy3.energy - 1.0/280.0*trial_energy4.energy)/delta_r;
                pair_force_LJ = -dFi_dr/r_0;
                fileOutput1 << r_0*(1e+10) << " " << molA.phi << " " << molB.phi << " " << pair_force_LJ <<  endl;

                dFi_dr = (1.0/280.0*trial_energy_4.energy_QQ - 4.0/105.0*trial_energy_3.energy_QQ + 1.0/5.0*trial_energy_2.energy_QQ - 4.0/5.0*trial_energy_1.energy_QQ + 4.0/5.0*trial_energy1.energy_QQ - 1.0/5.0*trial_energy2.energy_QQ + 4.0/105.0*trial_energy3.energy_QQ - 1.0/280.0*trial_energy4.energy_QQ)/delta_r;
                pair_force_QQ = -dFi_dr/r_0;
                fileOutput2 << r_0*(1e+10) << " " << molA.phi << " " << molB.phi << " " << pair_force_QQ <<  endl;
            }
        }

    }
     fileOutput.close();
     fileOutput1.close();
     fileOutput2.close();
}
