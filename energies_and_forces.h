results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    double distance = 0;

    results en_and_press;

    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
    double x1;
    for (int id = -1; id < 2; id++)
    {
       x1 = (molB_x - molA_x + id*Lx);
       if (abs(x1) > Rc) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          y1 = (molB_y - molA_y + jd*Ly);
          if (abs(y1) > Rc) {continue;}
             r2 = x1*x1 + y1*y1;
             if (r2 <= Rc2)
             {
               //calculation of energy for two molecules
             }
       }
    }

    energy *= beta;

    en_and_press.energy = energy;
    return en_and_press;
}
