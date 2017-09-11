void inter_for_pressure_ab(state molA, state molB, double &Rc, double &Rc2, double &Lx, double &Ly, const double &A, double &C_q,
                         double &beta, double &p_N, double &p_T, double p_Tot)
{
    //double sigma = 331.8e-12;                 // Sigma in pm
    //double dn2 = 109.8e-12;                   // Distance between nitrogen atoms in pm
    //double dq1 = 84.7e-12;                    // Distance between "+" charge and center of quadrupole in pm
    //double dq2 = 104.4e-12;                   // Distance between "-" charge and center of quadrupole in pm

    double dn2 = 0.33092224232;               // Distance between nitrogen atoms in sigma units
    double dq1 = 0.25527426160;               // Distance between "+" charge and center of quadrupole in sigma units
    double dq2 = 0.31464737794;               // Distance between "-" charge and center of quadrupole in sigma units

    double eps = 0.515e-21;                         // LJ energy for nitrogen in J
    const double qe = 1.6021766208e-19;                   // The charge of one electron in C
    double q = 0.373*qe;                            // Charge of the quadrupole points in C
    double q2 = q*q;
    double sigma = 0.3318e-9;
    double dist,dist2,a,b,c,vir_LJ,vir_QQ;

    valarray<double> l_i={sin(molA.tetta)*cos(molA.phi), sin(molA.tetta)*sin(molA.phi), cos(molA.tetta)};
    valarray<double> l_j={sin(molB.tetta)*cos(molB.phi), sin(molB.tetta)*sin(molB.phi), cos(molB.tetta)};
    valarray<double> r_ij;//={distPBC(Lx,molB.x - molA.x), distPBC(Ly,molB.y - molA.y), 0};

    double P_LJ_N=0, P_LJ_appr_N=0;
    double P_QQ_N=0, P_QQ_appr_N=0;
    double P_LJ_T=0, P_LJ_appr_T=0;
    double P_QQ_T=0, P_QQ_appr_T=0;
    double P_LJ=0, P_LJ_appr=0;
    double P_QQ=0, P_QQ_appr=0;
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
                r_ij = {x1, y1, 0};
                x1 = abs(x1);
                y1 = abs(y1);
        //////////////////////////////////////////////////////////
        ////////CALCULATION OF LJ VIRIAL PRESSURE FOR N2
        //////////////////////////////////////////////////////////

                // Exact calculation of LJ pressure

                //AC
                vect = dn2/2*l_i+r_ij-dn2/2*l_j;
                dist2 = (vect*vect).sum();
                double invDr6 = 1.0/pow(dist2, 3);
                dist = sqrt((vect*vect).sum());
                vir_LJ = invDr6 * (2*invDr6 - 1)/dist;
                P_LJ_N += vir_LJ*(abs(vect[0])/dist)*x1;
                P_LJ_T += vir_LJ*(abs(vect[1])/dist)*y1;
                P_LJ += vir_LJ*dist2;
                //BD
                vect = -dn2/2*l_i+r_ij+dn2/2*l_j;
                dist2 = (vect*vect).sum();
                invDr6 = 1.0/pow(dist2, 3);
                dist = sqrt((vect*vect).sum());
                vir_LJ = invDr6 * (2*invDr6 - 1)/dist;
                P_LJ_N += vir_LJ*(abs(vect[0])/dist)*x1;
                P_LJ_T += vir_LJ*(abs(vect[1])/dist)*y1;
                P_LJ += vir_LJ*dist2;
                //AD
                vect = dn2/2*l_i+r_ij+dn2/2*l_j;
                dist2 = (vect*vect).sum();
                invDr6 = 1.0/pow(dist2, 3);
                dist = sqrt((vect*vect).sum());
                vir_LJ = invDr6 * (2*invDr6 - 1)/dist;
                P_LJ_N += vir_LJ*(abs(vect[0])/dist)*x1;
                P_LJ_T += vir_LJ*(abs(vect[1])/dist)*y1;
                P_LJ += vir_LJ*dist2;
                //BC
                vect = -dn2/2*l_i+r_ij-dn2/2*l_j;
                dist2 = (vect*vect).sum();
                invDr6 = 1.0/pow(dist2, 3);
                dist = sqrt((vect*vect).sum());
                vir_LJ = invDr6 * (2*invDr6 - 1)/dist;
                P_LJ_N += vir_LJ*(abs(vect[0])/dist)*x1;
                P_LJ_T += vir_LJ*(abs(vect[1])/dist)*y1;
                P_LJ += vir_LJ*dist2;

                // Approximate calculation of the LJ interaction
                // through the point potential
            /*
                a=(r_ij*l_i).sum();
                b=(r_ij*l_j).sum();
                c=(l_i*l_j).sum();
                dist2=(r_ij*r_ij).sum();

                double s[4];
                s[0]=-3*dn2*(a-b)/dist2+3/2*pow(dn2,2)*(pow((a-b),2)/pow(dist2,2)-(1-c)/dist2);
                s[1]=3*dn2*(a-b)/dist2+3/2*pow(dn2,2)*(pow((a-b),2)/pow(dist2,2)-(1-c)/dist2);
                s[2]=-3*dn2*(a+b)/dist2+3/2*pow(dn2,2)*(pow((a+b),2)/pow(dist2,2)-(1+c)/dist2);
                s[3]=3*dn2*(a+b)/dist2+3/2*pow(dn2,2)*(pow((a+b),2)/pow(dist2,2)-(1+c)/dist2);

                double g1=0,g2=0,temp_g1;
                for (int i=0;i<4;i++) {temp_g1=exp(s[i]);g1+=temp_g1;g2+=pow(temp_g1,2);}

                U_LJ+=g2/pow(dist2,6)-g1/pow(dist2,3);
            */

            ///////////////////////////////////////////////////////////////
            ////////CALCULATION OF QQ INTERACTION OF TWO LINEAR QUADRUPOLES
            ///////////////////////////////////////////////////////////////

                // Exact calculation of QQ interaction
                // in A1B1C1D1 - A2B2C2D2 pair

                // A1A2
                vect = dq2*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // A1B2
                vect = dq2*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // A1C2
                vect = dq2*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // A1D2
                vect = dq2*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // B1A2
                vect = dq1*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // B1B2
                vect = dq1*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // B1C2
                vect = dq1*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // B1D2
                vect = dq1*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // C1A2
                vect = -dq1*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // C1B2
                vect = -dq1*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // C1C2
                vect = -dq1*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // C1D2
                vect = -dq1*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // D1A2
                vect = -dq2*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // D1B2
                vect = -dq2*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // D1C2
                vect = -dq2*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N -= vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T -= vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ -= vir_QQ*dist2;

                // D1D2
                vect = -dq2*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                vir_QQ = A*q2/dist*log(dist);
                P_QQ_N += vir_QQ*(abs(vect[0])/dist)*x1;
                P_QQ_T += vir_QQ*(abs(vect[1])/dist)*y1;
                P_QQ += vir_QQ*dist2;

                // Exact calculation of QQ interaction
                // in A1B1C1D1 - A2B2C2D2 pair
            /*
                a=(r_ij*l_i).sum();
                b=(r_ij*l_j).sum();
                c=(l_i*l_j).sum();
                dist2=(r_ij*r_ij).sum();
                U_QQ_appr += C_q*(1+2*pow(c,2)-5*(pow(a,2)+pow(b,2)+4*a*b*c)/dist2+35*pow(a*b,2)/pow(dist2,2))/pow(dist2,2.5);
            */
            //cout << "r_ij: " << sqrt(pow(r_ij[0],2)+pow(r_ij[1],2)+pow(r_ij[2],2))/sigma /*<< " LJ: " << U_LJ << " "*/ << "LJ: " << U_LJ << " QQ:" << U_QQ/eps << endl;


             }
       }
    }

    P_LJ_N *= 24*beta/Lx/Ly;
    P_LJ_T *= 24*beta/Lx/Ly;
    P_LJ *= 24*beta/Lx/Ly;
    P_QQ_N *= -beta/Lx/Ly/eps;
    P_QQ_T *= -beta/Lx/Ly/eps;
    P_QQ *= -beta/Lx/Ly/eps;


   p_N += P_LJ_N+P_QQ_N;
   p_T += P_LJ_T+P_QQ_T;
   p_Tot += P_LJ + P_QQ;
}
