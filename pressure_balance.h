void pressure_balance(double press_X, double press_Y, double &Lx, double &Ly, int &nPart, vector <state> &coordinates, double &beta)
{
    double dL,Lx_new,Ly_new;    // Linear sizes of the box after correction
    double S = Lx*Ly;
    dL = 0.005;
    double pressure_dif = abs((press_X-press_Y)/Lx/Ly/sigma/sigma*1000);

    if (pressure_dif >= 2.0) {BALANCE_STEPS = 100;}
    if (pressure_dif < 2.0 && pressure_dif >= 1.5) {BALANCE_STEPS = 200;}
    if (pressure_dif < 1.5 && pressure_dif >= 1.25) {BALANCE_STEPS = 300;}
    if (pressure_dif < 1.25 && pressure_dif >= 1.0) {BALANCE_STEPS = 500;}
    if (pressure_dif < 1.0 && pressure_dif >= 0.7) {BALANCE_STEPS = 1000;}
    if (pressure_dif < 0.7 && pressure_dif >= 0.5) {BALANCE_STEPS = 1500;}
    if (pressure_dif < 0.5) {BALANCE_STEPS = 2500;}

    cout << "press_dif=" << pressure_dif << " balance_steps=" << BALANCE_STEPS << endl;

       if(press_X < press_Y)
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
       PotentialEnergy(nPart, Lx, Ly, coordinates, beta);
       double AR = ACCEPTANCE_RATIO[1]/(ACCEPTANCE_RATIO[0]+ACCEPTANCE_RATIO[1]);
       if (AR < 0.25 && delta_angle > 10.0*(3.141592653589/180.0) && delta_angle <= 90.0*(3.141592653589/180.0))
        {delta_angle -= 5.0*(3.141592653589/180.0);}
       if (AR > 0.3 && delta_angle >= 5.0*(3.141592653589/180.0) && delta_angle < 85.0*(3.141592653589/180.0))
        {delta_angle += 5.0*(3.141592653589/180.0);}
}
