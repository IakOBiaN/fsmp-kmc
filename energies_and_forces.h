results energies_and_forces(state molA, state molB, double &Lx, double &Ly, double &beta)
{
    double q2 = q*q;
    double dist,dist2,b1,b2,g,invDr6,h1,h2,h3,vir_LJ,vir_QQ;

    valarray<double> l_i={cos(molA.phi), sin(molA.phi)};
    valarray<double> l_j={cos(molB.phi), sin(molB.phi)};
    valarray<double> r_ij;

    double U_LJ=0;
    double U_QQ=0;
    results en_and_press;
    en_and_press.p.X_LJ = 0;
    en_and_press.p.X_QQ = 0;
    en_and_press.p.Y_LJ = 0;
    en_and_press.p.Y_QQ = 0;

    valarray<double> vect;

    double xxx  = molA.x;
    double yyy  = molA.y;
    double x1, y1, r2;
    double x2 = molB.x;
    double y2 = molB.y;
    for (int id = -1; id < 2; id++)
    {
       x1 = (x2 - xxx + id*Lx);
       if (abs(x1)>Rc) {continue;}
       for (int jd = -1; jd < 2; jd++)
       {
          y1 = (y2 - yyy + jd*Ly);
          if (abs(y1)>Rc) {continue;}
             r2 = x1*x1 + y1*y1;
             if (r2 <= Rc2)
             {
                        r_ij = {x1, y1};

                        //////////////////////////////////////////////////////////
                        ////////CALCULATION OF LJ INTERACTION OF DIATOMIC MOLECULE
                        //////////////////////////////////////////////////////////

                        //Exact calculation of LJ interaction and/or pressure in AB - CD pair
                        //AC
                        vect = dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                        dist2 = (vect*vect).sum();
                        invDr6 = 1.0/pow(dist2, 3);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //BD
                        vect = -dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                        dist2 = (vect*vect).sum();
                        invDr6 = 1.0/pow(dist2, 3);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //AD
                        vect = dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                        dist2 = (vect*vect).sum();
                        invDr6 = 1.0/pow(dist2, 3);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //BC
                        vect = -dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                        dist2 = (vect*vect).sum();
                        invDr6 = 1.0/pow(dist2, 3);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;

                ////////////////////////////////////////////////////////////////
                //////// CALCULATION OF QQ INTERACTION OF TWO LINEAR QUADRUPOLES
                ////////////////////////////////////////////////////////////////

                if (energy_QQ_exact || pressure_QQ_exact)
                {
                        //Exact calculation of QQ interaction in A1B1C1D1 - A2B2C2D2 pair
                        // A1A2
                        vect = dq2*l_i+r_ij-dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // A1B2
                        vect = dq2*l_i+r_ij-dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // A1C2
                        vect = dq2*l_i+r_ij+dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // A1D2
                        vect = dq2*l_i+r_ij+dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1A2
                        vect = dq1*l_i+r_ij-dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // B1B2
                        vect = dq1*l_i+r_ij-dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1C2
                        vect = dq1*l_i+r_ij+dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1D2
                        vect = dq1*l_i+r_ij+dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // C1A2
                        vect = -dq1*l_i+r_ij-dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // C1B2
                        vect = -dq1*l_i+r_ij-dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // C1C2
                        vect = -dq1*l_i+r_ij+dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // C1D2
                        vect = -dq1*l_i+r_ij+dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1A2
                        vect = -dq2*l_i+r_ij-dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // D1B2
                        vect = -dq2*l_i+r_ij-dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1C2
                        vect = -dq2*l_i+r_ij+dq1*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1D2
                        vect = -dq2*l_i+r_ij+dq2*l_j;
                        dist2 = (vect*vect).sum();
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }
                }
                else
                {
                        //Approximate calculation of QQ interaction and pressure in A1B1C1D1 - A2B2C2D2 pair
                        b1 = (r_ij*l_i).sum();
                        b2 = (r_ij*l_j).sum();
                        g = (l_i*l_j).sum();
                        dist2 = (r_ij*r_ij).sum();

                        if (!energy_QQ_exact)
                        {
                            U_QQ += C_q*(1.0+2.0*pow(g,2)-5.0*(pow(b1,2)+pow(b2,2)+4.0*b1*b2*g)/dist2+35.0*pow((b1*b2),2)/pow(dist2,2))/pow(dist2,2.5);
                        }

                        if (!pressure_QQ_exact)
                        {
                            dist = sqrt(dist2);
                            h1 = 1.0 + 2.0*pow(g,2) - 7.0*(pow(b1,2) + pow(b2,2) + 4.0*g*b1*b2)/dist2 + 63.0*pow(b1,2)*pow(b2,2)/pow(dist2,2);
                            h2 = b1 + 2.0*g*b2 - 7.0*b1*pow(b2,2)/dist2;
                            h3 = b2 + 2.0*g*b1 - 7.0*pow(b1,2)*b2/dist2;
                            en_and_press.p.X_QQ += 5*C_q*r_ij[0]/pow(dist2,3)/dist*(h1*r_ij[0] + 2.0*(h2*l_i[0] + h3*l_j[0]));
                            en_and_press.p.Y_QQ += 5*C_q*r_ij[1]/pow(dist2,3)/dist*(h1*r_ij[1] + 2.0*(h2*l_i[1] + h3*l_j[1]));
                        }
                }

             }
       }
    }

    U_LJ *= 4.0*beta;
    U_QQ *= beta*sigma/eps;
    en_and_press.p.X_LJ *= 24.0*eps;
    en_and_press.p.Y_LJ *= 24.0*eps;
    en_and_press.p.X_QQ *= sigma;
    en_and_press.p.Y_QQ *= sigma;

    en_and_press.energy = U_LJ + U_QQ;
    //if(energy > gm){energy = gm;}
    return en_and_press;
}
