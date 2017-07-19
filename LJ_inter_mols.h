double LJ_inter_mols(state molA, state molB, double &Lx, double &Ly)
{
    double dn2 = 0.1098;                            // Distance between nitrogen atoms in nm
    //double Qn2 = -4.453e-40;                        // Quadrupole moment of N2 molecule
    double dq1 = 84.7e-12;                          // Distance between "+" charge and center of quadrupole in pm
    double dq2 = 104.4e-12;                         // Distance between "-" charge and center of quadrupole in pm
    double eps = 0.515e-21;                         // LJ energy for nitrogen in J
    double qe = 1.6021766208e-19;                   // The charge of one electron in C
    double q = 0.373*qe;                            // Charge of the quadrupole points in C
    double q2 = q*q;
    double eps0 = 8.85418781762e-12;                // The permittivity of free space in C2 m-2 N-1
    const double A = 1.0/(4.0*3.1415926535*eps0)/4.0;   // Coulomb's constant, *4 to meet LJ calculation
    double sigma = 0.3318e-9;

    double l_i[3]={sin(molA.tetta)*cos(molA.phi), sin(molA.tetta)*sin(molA.phi), cos(molA.tetta)};
    double l_j[3]={sin(molB.tetta)*cos(molB.phi), sin(molB.tetta)*sin(molB.phi), cos(molB.tetta)};
    double r_ij[3] = {distPBC(Lx,molB.x - molA.x), distPBC(Ly,molB.y - molA.y), 0};

    double U_LJ;

    //AC
    double vect[3] = {((dn2/2)*l_i[0]+r_ij[0]-(dn2/2)*l_j[0]),((dn2/2)*l_i[1]+r_ij[1]-(dn2/2)*l_j[1]),((dn2/2)*l_i[2]+r_ij[2]-(dn2/2)*l_j[2])};
    double dist2 = pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2);
    double invDr6 = 1.0/pow(dist2, 3);
    U_LJ = (invDr6 * (invDr6 - 1));
    //BD
    vect[0] = (-(dn2/2)*l_i[0]+r_ij[0]+(dn2/2)*l_j[0]);
    vect[1] = (-(dn2/2)*l_i[1]+r_ij[1]+(dn2/2)*l_j[1]);
    vect[2] = (-(dn2/2)*l_i[2]+r_ij[2]+(dn2/2)*l_j[2]);
    dist2 = pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2);
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));
    //AD
    vect[0] = ((dn2/2)*l_i[0]+r_ij[0]+(dn2/2)*l_j[0]);
    vect[1] = ((dn2/2)*l_i[1]+r_ij[1]+(dn2/2)*l_j[1]);
    vect[2] = ((dn2/2)*l_i[2]+r_ij[2]+(dn2/2)*l_j[2]);
    dist2 = pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2);
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));
    //BC
    vect[0] = (-(dn2/2)*l_i[0]+r_ij[0]-(dn2/2)*l_j[0]);
    vect[1] = (-(dn2/2)*l_i[1]+r_ij[1]-(dn2/2)*l_j[1]);
    vect[2] = (-(dn2/2)*l_i[2]+r_ij[2]-(dn2/2)*l_j[2]);
    dist2 = pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2);
    invDr6 = 1.0/pow(dist2, 3);
    U_LJ += (invDr6 * (invDr6 - 1));

    double U_QQ;
    double dist;
    r_ij[0] = distPBC(Lx,molB.x - molA.x)*sigma;
    r_ij[1] = distPBC(Ly,molB.y - molA.y)*sigma;
    r_ij[2] = 0;

    // A1A2
    vect[0] = dq2*l_i[0]+r_ij[0]-dq2*l_j[0];
    vect[1] = dq2*l_i[1]+r_ij[1]-dq2*l_j[1];
    vect[2] = dq2*l_i[2]+r_ij[2]-dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ = A*q2/dist;

    // A1B2
    vect[0] = dq2*l_i[0]+r_ij[0]-dq1*l_j[0];
    vect[1] = dq2*l_i[1]+r_ij[1]-dq1*l_j[1];
    vect[2] = dq2*l_i[2]+r_ij[2]-dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // A1C2
    vect[0] = dq2*l_i[0]+r_ij[0]+dq1*l_j[0];
    vect[1] = dq2*l_i[1]+r_ij[1]+dq1*l_j[1];
    vect[2] = dq2*l_i[2]+r_ij[2]+dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // A1D2
    vect[0] = dq2*l_i[0]+r_ij[0]+dq2*l_j[0];
    vect[1] = dq2*l_i[1]+r_ij[1]+dq2*l_j[1];
    vect[2] = dq2*l_i[2]+r_ij[2]+dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // B1A2
    vect[0] = dq1*l_i[0]+r_ij[0]-dq2*l_j[0];
    vect[1] = dq1*l_i[1]+r_ij[1]-dq2*l_j[1];
    vect[2] = dq1*l_i[2]+r_ij[2]-dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // B1B2
    vect[0] = dq1*l_i[0]+r_ij[0]-dq1*l_j[0];
    vect[1] = dq1*l_i[1]+r_ij[1]-dq1*l_j[1];
    vect[2] = dq1*l_i[2]+r_ij[2]-dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // B1C2
    vect[0] = dq1*l_i[0]+r_ij[0]+dq1*l_j[0];
    vect[1] = dq1*l_i[1]+r_ij[1]+dq1*l_j[1];
    vect[2] = dq1*l_i[2]+r_ij[2]+dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // B1D2
    vect[0] = dq1*l_i[0]+r_ij[0]+dq2*l_j[0];
    vect[1] = dq1*l_i[1]+r_ij[1]+dq2*l_j[1];
    vect[2] = dq1*l_i[2]+r_ij[2]+dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // C1A2
    vect[0] = -dq1*l_i[0]+r_ij[0]-dq2*l_j[0];
    vect[1] = -dq1*l_i[1]+r_ij[1]-dq2*l_j[1];
    vect[2] = -dq1*l_i[2]+r_ij[2]-dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // C1B2
    vect[0] = -dq1*l_i[0]+r_ij[0]-dq1*l_j[0];
    vect[1] = -dq1*l_i[1]+r_ij[1]-dq1*l_j[1];
    vect[2] = -dq1*l_i[2]+r_ij[2]-dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // C1C2
    vect[0] = -dq1*l_i[0]+r_ij[0]+dq1*l_j[0];
    vect[1] = -dq1*l_i[1]+r_ij[1]+dq1*l_j[1];
    vect[2] = -dq1*l_i[2]+r_ij[2]+dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // C1D2
    vect[0] = -dq1*l_i[0]+r_ij[0]+dq2*l_j[0];
    vect[1] = -dq1*l_i[1]+r_ij[1]+dq2*l_j[1];
    vect[2] = -dq1*l_i[2]+r_ij[2]+dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // D1A2
    vect[0] = -dq2*l_i[0]+r_ij[0]-dq2*l_j[0];
    vect[1] = -dq2*l_i[1]+r_ij[1]-dq2*l_j[1];
    vect[2] = -dq2*l_i[2]+r_ij[2]-dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    // D1B2
    vect[0] = -dq2*l_i[0]+r_ij[0]-dq1*l_j[0];
    vect[1] = -dq2*l_i[1]+r_ij[1]-dq1*l_j[1];
    vect[2] = -dq2*l_i[2]+r_ij[2]-dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // D1C2
    vect[0] = -dq2*l_i[0]+r_ij[0]+dq1*l_j[0];
    vect[1] = -dq2*l_i[1]+r_ij[1]+dq1*l_j[1];
    vect[2] = -dq2*l_i[2]+r_ij[2]+dq1*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ -= A*q2/dist;

    // D1D2
    vect[0] = -dq2*l_i[0]+r_ij[0]+dq2*l_j[0];
    vect[1] = -dq2*l_i[1]+r_ij[1]+dq2*l_j[1];
    vect[2] = -dq2*l_i[2]+r_ij[2]+dq2*l_j[2];
    dist = sqrt(pow(vect[0],2)+pow(vect[1],2)+pow(vect[2],2));
    U_QQ += A*q2/dist;

    //cout << "r_ij: " << sqrt(pow(r_ij[0],2)+pow(r_ij[1],2)+pow(r_ij[2],2))/sigma << " " << "LJ: " << U_LJ << " " << "QQ: " << U_QQ/eps << endl;

    return (U_LJ + U_QQ/eps);
}
