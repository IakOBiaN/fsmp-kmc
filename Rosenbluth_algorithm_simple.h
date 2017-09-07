// Choose a particle to be displaced
// according to the ROSENBLUTH scheme
int Rosenbluth_algorithm_simple(int &nPart, vector <state> &coordinates, double &dt)
{
 int trialPart = 0;

 coordinates[0].mob = exp(coordinates[0].energy);
 for(int i = 1; i < nPart; i++){
                                    coordinates[i].mob = coordinates[i-1].mob + exp(coordinates[i].energy);
                               }

 dt = 1.0/coordinates[nPart-1].mob;

 double Rp = coordinates[nPart-1].mob*RanGen.Random();

 if(Rp < coordinates[0].mob){trialPart = 0;}
   else
        {for(int i = 1; i < nPart; i++){if(Rp >= coordinates[i-1].mob && Rp < coordinates[i].mob){trialPart = i; break;}}}

 return trialPart;
}
