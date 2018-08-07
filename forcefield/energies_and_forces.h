results energies_and_forces(state molA, state molB, double Lx, double Ly)
{
    bool rosenbluth = false; //kMC NOT WORKING NOW!!!// If rosenbluth = false then Metropolis algorithm works
    bool energy_QQ_exact = true;
    bool pressure_QQ_exact = true;

    double q2 = q*q;
    double dist,dist2,dist4,b1,b1_2,b2_2,b2,b1b2,g,g_2,invDr6,h1,h2,h3,vir_LJ,vir_QQ;

    double l_i[2]={molA.cos_phi, molA.sin_phi};
    double l_j[2]={molB.cos_phi, molB.sin_phi};
    double r_ij[2], vect[2];

    //cout << l_i[0] << " " << l_i[1] << endl;

    double U_LJ=0;
    double U_QQ=0;
    results en_and_press;

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

                        r_ij[0] = x1;
                        r_ij[1] = y1;

                        //////////////////////////////////////////////////////////
                        ////////CALCULATION OF LJ INTERACTION OF DIATOMIC MOLECULE
                        //////////////////////////////////////////////////////////

                        //Exact calculation of LJ interaction and/or pressure in AB - CD pair
                        //AC
                        vect[0] = dn2/2.0*l_i[0]+r_ij[0]-dn2/2.0*l_j[0];
                        vect[1] = dn2/2.0*l_i[1]+r_ij[1]-dn2/2.0*l_j[1];
                        dist2 = (vect[0]*vect[0]+vect[1]*vect[1])/pow(sigma, 2);
                        invDr6 = 1.0/(dist2*dist2*dist2);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //BD
                        vect[0] = -dn2/2.0*l_i[0]+r_ij[0]+dn2/2.0*l_j[0];
                        vect[1] = -dn2/2.0*l_i[1]+r_ij[1]+dn2/2.0*l_j[1];
                        dist2 = (vect[0]*vect[0]+vect[1]*vect[1])/pow(sigma, 2);
                        invDr6 = 1.0/(dist2*dist2*dist2);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //AD
                        vect[0] = dn2/2.0*l_i[0]+r_ij[0]+dn2/2.0*l_j[0];
                        vect[1] = dn2/2.0*l_i[1]+r_ij[1]+dn2/2.0*l_j[1];
                        dist2 = (vect[0]*vect[0]+vect[1]*vect[1])/pow(sigma, 2);
                        invDr6 = 1.0/(dist2*dist2*dist2);
                        U_LJ += (invDr6 * (invDr6 - 1.0));
                        vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                        en_and_press.p.X_LJ += vir_LJ*vect[0]*x1;
                        en_and_press.p.Y_LJ += vir_LJ*vect[1]*y1;
                        //BC
                        vect[0] = -dn2/2.0*l_i[0]+r_ij[0]-dn2/2.0*l_j[0];
                        vect[1] = -dn2/2.0*l_i[1]+r_ij[1]-dn2/2.0*l_j[1];
                        dist2 = (vect[0]*vect[0]+vect[1]*vect[1])/pow(sigma, 2);
                        invDr6 = 1.0/(dist2*dist2*dist2);
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
                        vect[0] = dq2*l_i[0] + r_ij[0] - dq2*l_j[0];
						vect[1] = dq2*l_i[1] + r_ij[1] - dq2*l_j[1];
                        dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // A1B2
						vect[0] = dq2*l_i[0] + r_ij[0] - dq1*l_j[0];
						vect[1] = dq2*l_i[1] + r_ij[1] - dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // A1C2
						vect[0] = dq2*l_i[0] + r_ij[0] + dq1*l_j[0];
						vect[1] = dq2*l_i[1] + r_ij[1] + dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // A1D2
						vect[0] = dq2*l_i[0] + r_ij[0] + dq2*l_j[0];
						vect[1] = dq2*l_i[1] + r_ij[1] + dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1A2
						vect[0] = dq1*l_i[0] + r_ij[0] - dq2*l_j[0];
						vect[1] = dq1*l_i[1] + r_ij[1] - dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // B1B2
						vect[0] = dq1*l_i[0] + r_ij[0] - dq1*l_j[0];
						vect[1] = dq1*l_i[1] + r_ij[1] - dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1C2
                        vect[0] = dq1*l_i[0]+r_ij[0]+dq1*l_j[0];
						vect[1] = dq1*l_i[1] + r_ij[1] + dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // B1D2
						vect[0] = dq1*l_i[0] + r_ij[0] + dq2*l_j[0];
						vect[1] = dq1*l_i[1] + r_ij[1] + dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // C1A2
						vect[0] = -dq1*l_i[0] + r_ij[0] - dq2*l_j[0];
						vect[1] = -dq1*l_i[1] + r_ij[1] - dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // C1B2
						vect[0] = -dq1*l_i[0] + r_ij[0] - dq1*l_j[0];
						vect[1] = -dq1*l_i[1] + r_ij[1] - dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // C1C2
						vect[0] = -dq1*l_i[0] + r_ij[0] + dq1*l_j[0];
						vect[1] = -dq1*l_i[1] + r_ij[1] + dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // C1D2
						vect[0] = -dq1*l_i[0] + r_ij[0] + dq2*l_j[0];
						vect[1] = -dq1*l_i[1] + r_ij[1] + dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1A2
						vect[0] = -dq2*l_i[0] + r_ij[0] - dq2*l_j[0];
						vect[1] = -dq2*l_i[1] + r_ij[1] - dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }

                        // D1B2
						vect[0] = -dq2*l_i[0] + r_ij[0] - dq1*l_j[0];
						vect[1] = -dq2*l_i[1] + r_ij[1] - dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1C2
						vect[0] = -dq2*l_i[0] + r_ij[0] + dq1*l_j[0];
						vect[1] = -dq2*l_i[1] + r_ij[1] + dq1*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ -= A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ -= vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ -= vir_QQ*vect[1]*y1/dist;
                        }

                        // D1D2
						vect[0] = -dq2*l_i[0] + r_ij[0] + dq2*l_j[0];
						vect[1] = -dq2*l_i[1] + r_ij[1] + dq2*l_j[1];
						dist2 = vect[0] * vect[0] + vect[1] * vect[1];
                        dist = sqrt(dist2);
                        if (energy_QQ_exact) {U_QQ += A*q2/dist;}
                        if (pressure_QQ_exact)
                        {
                            vir_QQ =  A*q2/dist2;
                            en_and_press.p.X_QQ += vir_QQ*vect[0]*x1/dist;
                            en_and_press.p.Y_QQ += vir_QQ*vect[1]*y1/dist;
                        }
                }
                //Approximate calculation of QQ interaction and pressure in A1B1C1D1 - A2B2C2D2 pair

                        b1 = r_ij[0] * l_i[0] + r_ij[1] * l_i[1];
                        b2 = r_ij[0] * l_j[0] + r_ij[1] * l_j[1];
                        g = l_i[0] * l_j[0] + l_i[1] * l_j[1];
                        dist2 = r2;
                        dist = sqrt(dist2);
                        dist4 = dist2*dist2;
                        g_2 = g*g;
                        b1_2 = b1*b1;
                        b2_2 = b2*b2;
                        b1b2 = b1*b2;
                        if (!energy_QQ_exact)
                        {
                            U_QQ += C_q*(1.0 + 2.0*g_2 - 5.0*(b1_2 + b2_2 + 4.0*b1b2*g) / dist2 + 35.0*b1b2*b1b2/ dist4) / (dist4*dist);
                        }

                        if (!pressure_QQ_exact)
                        {
							h1 = 1.0 + 2.0*g_2 - 7.0*(b1_2 + b2_2 + 4.0*b1b2*g) / dist2 + 63.0*b1_2*b2_2 / dist4;							h2 = b1 + 2.0*g*b2 - 7.0*b1*b2_2 / dist2;							h3 = b2 + 2.0*g*b1 - 7.0*b1_2*b2 / dist2;							en_and_press.p.X_QQ += 5 * C_q*r_ij[0] / (dist4*dist2) / dist*(h1*r_ij[0] + 2.0*(h2*l_i[0] + h3*l_j[0]));							en_and_press.p.Y_QQ += 5 * C_q*r_ij[1] / (dist4*dist2) / dist*(h1*r_ij[1] + 2.0*(h2*l_i[1] + h3*l_j[1]));
                        }


             }
       }
    }

    U_LJ *= 4.0*eps;
    U_QQ *= 1.0; //sigma;
    en_and_press.p.X_LJ *= 24.0*eps;
    en_and_press.p.Y_LJ *= 24.0*eps;
    en_and_press.p.X_QQ *= 1.0;//sigma;
    en_and_press.p.Y_QQ *= 1.0;//sigma;

    //cout << "Px_LJ: " << en_and_press.p.X_LJ << "  Py_LJ: " << en_and_press.p.Y_LJ << "  Px_QQ: " << en_and_press.p.X_QQ << "  Py_QQ: " << en_and_press.p.Y_QQ << endl;

    //cout << "Px: " << (en_and_press.p.X_LJ + en_and_press.p.X_QQ) << "  Py: " << (en_and_press.p.Y_LJ + en_and_press.p.Y_QQ) << endl;

    en_and_press.energy = U_LJ;
    en_and_press.energy_QQ = U_QQ;
    //if(energy > gm){energy = gm;}
    return en_and_press;
}
