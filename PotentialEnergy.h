using namespace std;

void PotentialEnergy(int &nPart, double &Lx, double &Ly, vector <state> &coordinates, double &beta)
{

 results en_and_press;
 results en_and_press_2;

EN_AND_PR_counter = en_and_press;

if (rosenbluth) {for(int i = 0; i < nPart; i++) {coordinates[i].energy = 0;}}

 // Loop over all distinct particle pairs
 for(int molA = 0; molA < (nPart - 1); molA++)
    {
     for(int molB = (molA + 1); molB < nPart; molB++)
      {
           en_and_press = energies_and_forces_2(coordinates[molA], coordinates[molB], Lx, Ly, beta);
           //en_and_press_2 = energies_and_forces(coordinates[molA], coordinates[molB], Lx, Ly, beta);

           /*if (abs(en_and_press.p_Y_LJ/en_and_press_2.p_Y_LJ) < 0.9 || abs(en_and_press.p_Y_LJ/en_and_press_2.p_Y_LJ) > 1.1)
           {
             cout << "exact_press=" << en_and_press.p_X_LJ << " " << en_and_press.p_X_QQ << " " << en_and_press.p_Y_LJ << " " << en_and_press.p_Y_QQ << endl;
             cout << "approx_press=" << en_and_press_2.p_X_LJ << " " << en_and_press_2.p_X_QQ << " " << en_and_press_2.p_Y_LJ << " " << en_and_press_2.p_Y_QQ << endl;
           }*/
           if (rosenbluth)
           {
                coordinates[molA].energy += en_and_press.energy;
                coordinates[molB].energy += en_and_press.energy;
           }
           EN_AND_PR_counter = EN_AND_PR_counter + en_and_press;
      }
    }
}
