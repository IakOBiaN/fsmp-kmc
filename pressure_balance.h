void pressure_balance(double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL, Lx_central_new, Ly_new;    // Linear sizes of the box after correction
    double Lx_central = Lx/8.0;
    double S = Lx_central*Ly; // Area of the central cell
    dL = dr/2.0; //0.005;
       if(press_X < press_Y)
       {
        Lx_central_new = Lx_central - dL;
        Ly_new = S/Lx_central_new;  // Area of the adlayer is kept constant
       }
       else
       {
        Lx_central_new = Lx_central + dL;
        Ly_new = S/Lx_central_new;
       }

       for(int i = 0; i < nPart; i++)
          {
            // Coordinates of the molecules change proportionally to
            // the change of the box size
            coordinates[i].x = (Lx_central_new/Lx_central)*coordinates[i].x;
            coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
          }
       Lx = Lx_central_new*8.0;
       Ly = Ly_new;
       // Recalculate energies after compressing or expanding the box
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
       double AR_r = ACCEPTANCE_RATIO_r[1]/(ACCEPTANCE_RATIO_r[0]+ACCEPTANCE_RATIO_r[1]);
       if (AR_r < 0.25 && delta_angle > 10.0)
        {delta_angle -= 5.0;}
       if (AR_r > 0.3 && delta_angle < 85.0)
        {delta_angle += 5.0;}
       double AR_m = ACCEPTANCE_RATIO_m[1]/(ACCEPTANCE_RATIO_m[0]+ACCEPTANCE_RATIO_m[1]);
       if (AR_m < 0.25 && delta > 0.05)
        {delta -= 0.05;}
       if (AR_m > 0.3 && delta < 0.8)
        {delta += 0.05;}
}
