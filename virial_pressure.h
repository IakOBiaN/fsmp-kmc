 void virial_pressure (int &nPart, double &Lx, double &Ly, double &beta, vector <state> &coordinates,
                         double &p_X_LJ, double &p_X_QQ, double &p_Y_LJ, double &p_Y_QQ)
 {

  double q2 = q*q;
  double dist,dist2,a,b,c,vir_LJ,vir_QQ;

  // Loop over all distinct particle pairs
  //double dr = 0.01;
  for (int i = 0; i < nPart-1; i++)
  {
     for (int j = i+1; j < nPart; j++)
     {
        //inter_for_pressure_ab(coordinates[molA], coordinates[molB], Rc, Rc2, Lx, Ly, A, C_q, beta, p_N, p_T, p_Tot);

        state molA = coordinates[i];
        state molB = coordinates[j];

        valarray<double> l_i={cos(molA.phi), sin(molA.phi)};
        valarray<double> l_j={cos(molB.phi), sin(molB.phi)};
        valarray<double> r_ij;

        double P_LJ_X=0;
        double P_QQ_X=0, P_QQ_appr_X=0;
        double P_LJ_Y=0;
        double P_QQ_Y=0, P_QQ_appr_Y=0;
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
            ////////CALCULATION OF LJ VIRIAL PRESSURE FOR N2
            //////////////////////////////////////////////////////////

                    // Exact calculation of LJ pressure

                    //AC
                    vect = dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                    dist2 = (vect*vect).sum();
                    double invDr6 = 1.0/pow(dist2, 3);
                    //if(dist2 < 0.56565441) {p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0; return;}
                    vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                    P_LJ_X += vir_LJ*vect[0]*x1;
                    P_LJ_Y += vir_LJ*vect[1]*y1;
                    //BD
                    vect = -dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                    dist2 = (vect*vect).sum();
                    invDr6 = 1.0/pow(dist2, 3);
                    //if(dist2 < 0.56565441) {p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0; return;}
                    vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                    P_LJ_X += vir_LJ*vect[0]*x1;
                    P_LJ_Y += vir_LJ*vect[1]*y1;
                    //AD
                    vect = dn2/2.0*l_i+r_ij+dn2/2.0*l_j;
                    dist2 = (vect*vect).sum();
                    invDr6 = 1.0/pow(dist2, 3);
                    //if(dist2 < 0.56565441) {p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0; return;}
                    vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                    P_LJ_X += vir_LJ*vect[0]*x1;
                    P_LJ_Y += vir_LJ*vect[1]*y1;
                    //BC
                    vect = -dn2/2.0*l_i+r_ij-dn2/2.0*l_j;
                    dist2 = (vect*vect).sum();
                    invDr6 = 1.0/pow(dist2, 3);
                    //if(dist2 < 0.56565441) {p_X_LJ = 0; p_X_QQ = 0; p_Y_LJ = 0; p_Y_QQ = 0; return;}
                    vir_LJ = invDr6 * (2.0*invDr6 - 1.0)/dist2;
                    P_LJ_X += vir_LJ*vect[0]*x1;
                    P_LJ_Y += vir_LJ*vect[1]*y1;

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
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // A1B2
                    vect = dq2*l_i+r_ij-dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // A1C2
                    vect = dq2*l_i+r_ij+dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // A1D2
                    vect = dq2*l_i+r_ij+dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // B1A2
                    vect = dq1*l_i+r_ij-dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // B1B2
                    vect = dq1*l_i+r_ij-dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // B1C2
                    vect = dq1*l_i+r_ij+dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // B1D2
                    vect = dq1*l_i+r_ij+dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // C1A2
                    vect = -dq1*l_i+r_ij-dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // C1B2
                    vect = -dq1*l_i+r_ij-dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // C1C2
                    vect = -dq1*l_i+r_ij+dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // C1D2
                    vect = -dq1*l_i+r_ij+dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // D1A2
                    vect = -dq2*l_i+r_ij-dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // D1B2
                    vect = -dq2*l_i+r_ij-dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // D1C2
                    vect = -dq2*l_i+r_ij+dq1*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X -= vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y -= vir_QQ*vect[1]*y1/dist;

                    // D1D2
                    vect = -dq2*l_i+r_ij+dq2*l_j;
                    dist2 = (vect*vect).sum();
                    dist = sqrt(dist2);
                    vir_QQ =  A*q2/dist2;
                    P_QQ_X += vir_QQ*vect[0]*x1/dist;
                    P_QQ_Y += vir_QQ*vect[1]*y1/dist;

                    // Approximate calculation of QQ interaction
                    // in A1B1C1D1 - A2B2C2D2 pair

                    double b1 = (r_ij*l_i).sum();
                    double b2 = (r_ij*l_j).sum();
                    double g = (l_i*l_j).sum();
                    dist2 = (r_ij*r_ij).sum();
                    dist = sqrt(dist2);
                    double h1 = 1.0 + 2.0*pow(g,2) - 7.0*(pow(b1,2) + pow(b2,2) + 4.0*g*b1*b2)/dist2 + 63.0*pow(b1,2)*pow(b2,2)/pow(dist2,2);
                    double h2 = b1 + 2.0*g*b2 - 7.0*b1*pow(b2,2)/dist2;
                    double h3 = b2 + 2.0*g*b1 - 7.0*pow(b1,2)*b2/dist2;

                    P_QQ_appr_X += 5*C_q*r_ij[0]/pow(dist2,3)/dist*(h1*r_ij[0] + 2.0*(h2*l_i[0] + h3*l_j[0]))*sigma;
                    P_QQ_appr_Y += 5*C_q*r_ij[1]/pow(dist2,3)/dist*(h1*r_ij[1] + 2.0*(h2*l_i[1] + h3*l_j[1]))*sigma;



                 }
           }
        }

        P_LJ_X *= 24.0*eps;
        P_LJ_Y *= 24.0*eps;
        //P_QQ_X *= 1.0*sigma;
        //P_QQ_Y *= 1.0*sigma;

        p_X_LJ += P_LJ_X;
        p_X_QQ += P_QQ_appr_X;
        p_Y_LJ += P_LJ_Y;
        p_Y_QQ += P_QQ_appr_Y;
     }
  }

}
