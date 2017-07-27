double Inter_potential(state molA, state molB, double &Lx, double &Ly, const double &A, double &C_q)
{
    double dn2 = 0.1098;                            // Distance between nitrogen atoms in nm
    double dq1 = 84.7e-12;                          // Distance between "+" charge and center of quadrupole in pm
    double dq2 = 104.4e-12;                         // Distance between "-" charge and center of quadrupole in pm
    double eps = 0.515e-21;                         // LJ energy for nitrogen in J
    const double qe = 1.6021766208e-19;                   // The charge of one electron in C
    double q = 0.373*qe;                            // Charge of the quadrupole points in C
    double q2 = q*q;
    double sigma = 0.3318e-9;
    double dist,dist2,a,b,c;

    valarray<double> l_i={sin(molA.tetta)*cos(molA.phi), sin(molA.tetta)*sin(molA.phi), cos(molA.tetta)};
    valarray<double> l_j={sin(molB.tetta)*cos(molB.phi), sin(molB.tetta)*sin(molB.phi), cos(molB.tetta)};
    valarray<double> r_ij={distPBC(Lx,molB.x - molA.x), distPBC(Ly,molB.y - molA.y), 0};

    double U_LJ, U_LJ_appr;
    double U_QQ, U_QQ_appr;

//////////////////////////////////////////////////////////
////////CALCULATION OF LJ INTERACTION OF DIATOMIC MOLECULE
//////////////////////////////////////////////////////////

    // Exact calculation of LJ interaction in AB - CD pair

    //AC
    valarray<double> vect = dn2/2*l_i+r_ij-dn2/2*l_j;
    dist2 = (vect*vect).sum();
    double invDr6 = 1.0/pow(dist2, 3);
    U_LJ = (invDr6 * (invDr6 - 1));
    //BD
    vect = -dn2/2*l_i+r_ij+dn2/2*l_j;
    dist2 = (vect*vect).sum();
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));
    //AD
    vect = dn2/2*l_i+r_ij+dn2/2*l_j;
    dist2 = (vect*vect).sum();
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));
    //BC
    vect = -dn2/2*l_i+r_ij-dn2/2*l_j;
    dist2 = (vect*vect).sum();
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));



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

    U_LJ=g2/pow(dist2,6)-g1/pow(dist2,3);
*/

///////////////////////////////////////////////////////////////
////////CALCULATION OF QQ INTERACTION OF TWO LINEAR QUADRUPOLES
///////////////////////////////////////////////////////////////

r_ij *= sigma; //correction for QQ interaction (distance in SI units)

    // Exact calculation of QQ interaction
    // in A1B1C1D1 - A2B2C2D2 pair

    // A1A2
    vect = dq2*l_i+r_ij-dq2*l_j;
    dist = sqrt((vect*vect).sum());
    U_QQ = A*q2/dist;

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

    // Exact calculation of QQ interaction
    // in A1B1C1D1 - A2B2C2D2 pair
/*
    a=(r_ij*l_i).sum();
    b=(r_ij*l_j).sum();
    c=(l_i*l_j).sum();
    dist2=(r_ij*r_ij).sum();
    U_QQ_appr = C_q*(1+2*pow(c,2)-5*(pow(a,2)+pow(b,2)+4*a*b*c)/dist2+35*pow(a*b,2)/pow(dist2,2))/pow(dist2,2.5);
*/
//cout << "r_ij: " << sqrt(pow(r_ij[0],2)+pow(r_ij[1],2)+pow(r_ij[2],2))/sigma /*<< " LJ: " << U_LJ << " "*/ << "LJ: " << U_LJ << " QQ:" << U_QQ/eps << endl;

    return (U_LJ + U_QQ/eps);
}
