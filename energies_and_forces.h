results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    results en_and_press;
    //cout << "START_EN_AND_FORCE" << endl;
    double ang1,ang2;
    double molA_x  = molA.x;
    double molA_y  = molA.y;
    double molB_x = molB.x;
    double molB_y = molB.y;
    double dx, dy,r2,energy = 0,pressure_X_LJ = 0, pressure_Y_LJ = 0, pressure_X_QQ = 0, pressure_Y_QQ = 0;
    double ang_molA = molA.phi;
    double ang_molB = molB.phi;
    bool mirror_int = true;
    ver v_000,v_010, v_100, v_110, v_001, v_011, v_101, v_111;
    coord v_find;
    int temp_dist,temp_ang1,temp_ang2;
    double dist_n;

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
               ang1 = ang_molA - dang;
               ang2 = ang_molB - dang;
               if (ang1<0) {ang1 += 360.0;}
               if (ang2<0) {ang2 += 360.0;}
               if (ang1>359.5) {ang1 -= 360.0;}
               if (ang2>359.5) {ang2 -= 360.0;}
               dist_n = (r-min_dist)/dr;
               dist = (int)(dist_n+0.5);
               a1 = (int)((ang1/da)+0.5);
               a2 = (int)((ang2/da)+0.5);
               if (r<min_dist)
               {
                 energy += forcefield[0][a1][a2]*100*exp(r/min_dist*log(0.01));
                 pressure_X_LJ += force_LJ[0][a1][a2]*100*exp(r/min_dist*log(0.01))*dx*dx;
                 pressure_Y_LJ += force_LJ[0][a1][a2]*100*exp(r/min_dist*log(0.01))*dy*dy;
                 pressure_X_QQ += force_QQ[0][a1][a2]*100*exp(r/min_dist*log(0.01))*dx*dx;
                 pressure_Y_QQ += force_QQ[0][a1][a2]*100*exp(r/min_dist*log(0.01))*dy*dy;
               }
               else
               {
                 //energy += forcefield[dist][a1][a2];
                 temp_dist = (int)(dist_n);
                 temp_ang1 = (int)(ang1/da);
                 temp_ang2 = (int)(ang2/da);
                 v_000.place.x = temp_dist;
                 v_000.place.y = temp_ang1;
                 v_000.place.z = temp_ang2;
                 v_000.value = forcefield[v_000.place.x][v_000.place.y][v_000.place.z];
                 v_010.place.x = temp_dist;
                 v_010.place.y = temp_ang1+1;
                 v_010.place.z = temp_ang2;
                 v_010.value = forcefield[v_010.place.x][v_010.place.y][v_010.place.z];
                 v_100.place.x = temp_dist+1;
                 v_100.place.y = temp_ang1;
                 v_100.place.z = temp_ang2;
                 v_100.value = forcefield[v_100.place.x][v_100.place.y][v_100.place.z];
                 v_110.place.x = temp_dist+1;
                 v_110.place.y = temp_ang1+1;
                 v_110.place.z = temp_ang2;
                 v_110.value = forcefield[v_110.place.x][v_110.place.y][v_110.place.z];
                 v_001.place.x = temp_dist;
                 v_001.place.y = temp_ang1;
                 v_001.place.z = temp_ang2+1;
                 v_001.value = forcefield[v_001.place.x][v_001.place.y][v_001.place.z];
                 v_011.place.x = temp_dist;
                 v_011.place.y = temp_ang1+1;
                 v_011.place.z = temp_ang2+1;
                 v_011.value = forcefield[v_011.place.x][v_011.place.y][v_011.place.z];
                 v_101.place.x = temp_dist+1;
                 v_101.place.y = temp_ang1;
                 v_101.place.z = temp_ang2+1;
                 v_101.value = forcefield[v_101.place.x][v_101.place.y][v_101.place.z];
                 v_111.place.x = temp_dist+1;
                 v_111.place.y = temp_ang1+1;
                 v_111.place.z = temp_ang2+1;
                 v_111.value = forcefield[v_111.place.x][v_111.place.y][v_111.place.z];
                 v_find.x = dist_n;
                 v_find.y = ang1/da;
                 v_find.z = ang2/da;
                 energy += interpol(v_000, v_010, v_100, v_110, v_001, v_011, v_101, v_111, v_find);
                 pressure_X_LJ += force_LJ[dist][a1][a2]*dx*dx;
                 pressure_Y_LJ += force_LJ[dist][a1][a2]*dy*dy;
                 pressure_X_QQ += force_QQ[dist][a1][a2]*dx*dx;
                 pressure_Y_QQ += force_QQ[dist][a1][a2]*dy*dy;
               }
             }
       }
    }
    en_and_press.energy = energy/temperature;
    en_and_press.p_X_LJ = pressure_X_LJ/temperature;
    en_and_press.p_Y_LJ = pressure_Y_LJ/temperature;
    en_and_press.p_X_QQ = pressure_X_QQ/temperature;
    en_and_press.p_Y_QQ = pressure_Y_QQ/temperature;
    return en_and_press;
}
