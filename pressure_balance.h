void pressure_balance(double &press_N, double &press_T, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &Rc,
                     double &Rc2, const double &A, double &C_q, double &beta)
{
    double dL,Lx_new,Ly_new;
    double S = Lx*Ly;
    dL = 0.005;

    if(abs(press_N - press_T) > 0.05)
    {
       if(press_T > press_N)
       {
        Lx_new = Lx - dL;
        Ly_new = S/Lx_new;
        for(int i = 0; i < nPart; i++)
           {
            coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
            coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
           }
        Lx = Lx_new;
        Ly = Ly_new;
       }
       else
       {
        Lx_new = Lx + dL;
        Ly_new = S/Lx_new;
        for(int i = 0; i < nPart; i++)
           {
            coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
            coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
           }
        Lx = Lx_new;
        Ly = Ly_new;
       }
       PotentialEnergy(nPart, Lx, Ly, Rc, Rc2, coordinates, A, C_q, beta);
    }
}
