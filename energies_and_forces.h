results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;

    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
    double dx, dy,r2,energy = 0;
    double ang1 = molA.phi;
    double ang2 = molB.phi;


    int dist;

    double r;
    for (int id = -1; id < 2; id++)
    {
       dx = (molB_x - molA_x + id*Lx);
       if (abs(dx) > max_dist) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          dy = (molB_y - molA_y + jd*Ly);
          if (abs(dy) > max_dist) {continue;}
             r2 = dx*dx + dy*dy;
             r = sqrt(r2);
             if (r < min_dist)
             {
               en_and_press.energy = 1e10;
               return en_and_press;
             }
             if (r <= max_dist)
             {
               //calculation of energy for two molecules
               double dang = acos(dx/r);
               ang1 = ang1 - dang;
               ang2 = ang2 - dang;
               if (ang1<0) {ang1 += 360.0;}
               if (ang2<0) {ang2 += 360.0;}
               if (ang1>360) {ang1 -= 360.0;}
               if (ang2>360) {ang2 -= 360.0;}
               dist = ceil((r-min_dist)/dr);
               //cout << "dist=" << dist << " ang1=" << ceil(ang1/da) << " ang2=" << ceil(ang2/da) << endl;
               energy += TMA_forcefield[dist][ceil(ang1/da)-1][ceil(ang2/da)-1];
             }
       }
    }

    energy *= beta;

    en_and_press.energy = energy;
    return en_and_press;
}
