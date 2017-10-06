void pressure_balance(double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &Rc,
                     double &Rc2, const double &A, double &C_q, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double S = Lx*Ly;
    dL = 0.005;
    double eps = 0.515e-21;

       if(press_X > press_Y)
       {
        Lx_new = Lx - dL;
        Ly_new = S/Lx_new;  // Area of the adlayer is kept constant
       }
       else
       {
        Lx_new = Lx + dL;
        Ly_new = S/Lx_new;
       }
       for(int i = 0; i < nPart; i++)
          {
            // Coordinates of the molecules change proportionally to
            // the change of the box size
            coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
            coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
          }
       Lx = Lx_new;
       Ly = Ly_new;
       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, Rc, Rc2, coordinates, A, C_q, beta);
}
