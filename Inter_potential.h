double Inter_potential(state molA, state molB, double &Rc, double &Rc2, double &Lx, double &Ly, const double &A, double &C_q, double &beta)
{
    //double sigma = 331.8e-12;                 // Sigma in pm
    //double dn2 = 109.8e-12;                   // Distance between nitrogen atoms in pm
    //double dq1 = 84.7e-12;                    // Distance between "+" charge and center of quadrupole in pm
    //double dq2 = 104.4e-12;                   // Distance between "-" charge and center of quadrupole in pm

    double dn2 = 0.33092224232;               // Distance between nitrogen atoms in sigma units
    double dq1 = 0.25527426160;               // Distance between "+" charge and center of quadrupole in sigma units
    double dq2 = 0.31464737794;               // Distance between "-" charge and center of quadrupole in sigma units


    double eps = 0.515e-21;                         // LJ energy for nitrogen in J
    const double qe = 1.6021766208e-19;             // The charge of one electron in C
    double q = 0.373*qe;                            // Charge of the quadrupole points in C
    double q2 = q*q;
    double dist,dist2,a,b,c;
    double gm = 50;

    valarray<double> l_i={cos(molA.phi), sin(molA.phi)};
    valarray<double> l_j={cos(molB.phi), sin(molB.phi)};
    valarray<double> r_ij;

    double U_LJ=0, U_LJ_appr=0;
    double U_QQ=0, U_QQ_appr=0;
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

                // Exact calculation of LJ interaction in AB - CD pair

                //AC
                vect = dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                dist2 = (vect*vect).sum();
                double invDr6 = 1.0/pow(dist2, 3);
                U_LJ += (invDr6 * (invDr6 - 1.0));
                //BD
                vect = -dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                dist2 = (vect*vect).sum();
                //cout << r_ij[1] << endl;
                invDr6 = 1.0/pow(dist2, 3);
                U_LJ += (invDr6 * (invDr6 - 1.0));
                //AD
                vect = dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                dist2 = (vect*vect).sum();
                invDr6 = 1.0/pow(dist2, 3);
                U_LJ += (invDr6 * (invDr6 - 1.0));
                //BC
                vect = -dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                dist2 = (vect*vect).sum();
                invDr6 = 1.0/pow(dist2, 3);
                U_LJ += (invDr6 * (invDr6 - 1.0));



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

            ////////////////////////////////////////////////////////////////
            //////// CALCULATION OF QQ INTERACTION OF TWO LINEAR QUADRUPOLES
            ////////////////////////////////////////////////////////////////

                // Exact calculation of QQ interaction
                // in A1B1C1D1 - A2B2C2D2 pair

                // A1A2
                vect = dq2*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // A1B2
                vect = dq2*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // A1C2
                vect = dq2*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // A1D2
                vect = dq2*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // B1A2
                vect = dq1*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // B1B2
                vect = dq1*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // B1C2
                vect = dq1*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // B1D2
                vect = dq1*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // C1A2
                vect = -dq1*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // C1B2
                vect = -dq1*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // C1C2
                vect = -dq1*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // C1D2
                vect = -dq1*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // D1A2
                vect = -dq2*l_i+r_ij-dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                // D1B2
                vect = -dq2*l_i+r_ij-dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // D1C2
                vect = -dq2*l_i+r_ij+dq1*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ -= A*q2/dist;

                // D1D2
                vect = -dq2*l_i+r_ij+dq2*l_j;
                dist = sqrt((vect*vect).sum());
                U_QQ += A*q2/dist;

                //cout << "r: " << dist << "\t" << "dn2: " << dn2 << endl;

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
    //cout << "LJ=" << 4*beta*U_LJ << " QQ=" << U_QQ/eps << endl;
    double energy = 4.0*beta*(U_LJ + 331.8e-12*U_QQ/eps/4.0);
    if(energy > gm){energy = gm;}
    return energy;
}
