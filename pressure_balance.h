void pressure_balance(double &press_N, double &press_T, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &Rc,
                     double &Rc2, const double &A, double &C_q, double &beta)
{
    double dL;
    if(abs(press_N - press_T) > 0.05)
    {
       dL = 0.005;
       if(press_T > press_N)
       {
        for(int i = 0; i < nPart; i++)
           {
            coordinates[i].x = ((Lx + dL)/Lx)*coordinates[i].x;
            coordinates[i].y = ((Ly + dL)/Ly)*coordinates[i].y;
           }
        Lx -= dL;
       }
       else
       {
        for(int i = 0; i < nPart; i++)
           {
            coordinates[i].x = ((Lx - dL)/Lx)*coordinates[i].x;
            coordinates[i].y = ((Ly - dL)/Ly)*coordinates[i].y;
           }
        Lx += dL;
       }
       PotentialEnergy(nPart, Lx, Ly, Rc, Rc2, coordinates, A, C_q, beta);
    }
}
