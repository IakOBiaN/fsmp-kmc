void layer_map(int &nPart, vector <state> &coordinates, vector <vector <double>> &xy_matrix, double &Lx, double &Ly)
{
    double dn2 = 0.5;

    for(int i = 0; i < nPart; i++)
    {
        double r = (1 - 2 * RanGen.Random());
        double dop = coordinates[i].x + dn2 * coordinates[i].cos_phi*r;
        if (dop > Lx) {dop -= Lx;}
        if (dop < 0)  {dop += Lx;}
        int xxx = dop*10;
        dop = coordinates[i].y + dn2 * coordinates[i].sin_phi*r;
        if (dop > Ly) {dop -= Ly;}
        if (dop < 0)  {dop += Ly;}
        int yyy = dop*10;
        xy_matrix[xxx][yyy]++;
    }
}
