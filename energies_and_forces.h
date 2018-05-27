results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
    //cout << "START_EN_AND_FORCE" << endl;
    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
    double dx, dy,r2,energy = 0;
    double ang1 = molA.phi;
    double ang2 = molB.phi;
    bool mirror_int = false;

    int from,to;
    if (mirror_int)
    {
    from = -1;
    to = 2;
    }
    else
    {
    from = 0;
    to = 1;
    }

    int dist,a1,a2;

    double r;
    for (int id = from; id < to; id++)
    {
       dx = (molB_x - molA_x + id*Lx);
       if (abs(dx) > max_dist) {continue;}
       for (int jd = from; jd < to; jd++)
       {
          dy = (molB_y - molA_y + jd*Ly);
          if (abs(dy) > max_dist) {continue;}
             r2 = dx*dx + dy*dy;
             r = sqrt(r2);
             if (r <= max_dist)
             {
               //calculation of energy for two molecules
               double dang = dx/r;
               if (dang<0) {dang=-acos(dang)/PI*180.0;} else {dang=acos(dang)/PI*180.0;}
               cout << "ACOS=" << dx/r << " and " << dang << endl;
               ang1 = ang1 - dang;
               ang2 = ang2 - dang;
               if (ang1<0) {ang1 += 360.0;}
               if (ang2<0) {ang2 += 360.0;}
               if (ang1>360) {ang1 -= 360.0;}
               if (ang2>360) {ang2 -= 360.0;}
               dist = (int)(((r-min_dist)/dr)+0.5);
               a1 = (int)((ang1/da)+0.5);
               a2 = (int)((ang2/da)+0.5);
               cout << "force=" << dist << " and " << ang1 << " and " << ang2 << endl;
               if (r<min_dist)
               {
                 energy += TMA_forcefield[0][a1][a2]*100*exp(r/min_dist*log(0.01));
               }
               else
               {
                 energy += TMA_forcefield[dist][a1][a2];
               }
             }
       }
    }
    en_and_press.energy = energy;
    return en_and_press;
}
