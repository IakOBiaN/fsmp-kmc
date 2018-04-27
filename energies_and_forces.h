results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    double distance = 0;

    results en_and_press;

    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
    double dx, dy;
    double ang1 = molA.ang;
    double ang2 = molB.ang;


    int dist;

    double r;
    for (int id = -1; id < 2; id++)
    {
       dx = (molB_x - molA_x + id*Lx);
       if (abs(dx) > Rc) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          dy = (molB_y - molA_y + jd*Ly);
          if (abs(dy) > Rc) {continue;}
             r2 = dx*dx + dy*dy;
             r = sqrt(r2);
             if (r <= Rc)
             {
               //calculation of energy for two molecules
               double dang = acos(dx/r);
               ang1 = ang1 - dang;
               ang2 = ang2 - dang;
               dist = ceil(r/dr);
               energy += forcefield[dist][ceil(ang1/da)][ceil(ang2/da)];
             }
       }
    }

    energy *= beta;

    en_and_press.energy = energy;
    return en_and_press;
}
