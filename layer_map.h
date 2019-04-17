void layer_map_TMA(int &nPart, vector <state> &coordinates, vector <vector <double> > &xy_matrix, double &Lx, double &Ly)
{
    double dn2 = 3.536;
    double r = RanGen.Random();
    double scale = 5.0;

    for(int i = 0; i < nPart; i++)
    {
        double dop = coordinates[i].x + dn2 * coordinates[i].cos_phi*r;
        if (dop > Lx) {dop -= Lx;}
        if (dop < 0)  {dop += Lx;}
        int xxx = dop*scale;
        dop = coordinates[i].y + dn2 * coordinates[i].sin_phi*r;
        if (dop > Ly) {dop -= Ly;}
        if (dop < 0)  {dop += Ly;}
        int yyy = dop*scale;
        xy_matrix[xxx][yyy]++;

        dop = coordinates[i].x + dn2 * (-0.5*coordinates[i].cos_phi - 0.86602540378443864676372317075294*coordinates[i].sin_phi)*r;
        if (dop > Lx) {dop -= Lx;}
        if (dop < 0)  {dop += Lx;}
        xxx = dop*scale;
        dop = coordinates[i].y + dn2 * (-0.5*coordinates[i].sin_phi + 0.86602540378443864676372317075294*coordinates[i].cos_phi)*r;
        if (dop > Ly) {dop -= Ly;}
        if (dop < 0)  {dop += Ly;}
        yyy = dop*scale;
        xy_matrix[xxx][yyy]++;

        dop = coordinates[i].x + dn2 * (-0.5*coordinates[i].cos_phi + 0.86602540378443864676372317075294*coordinates[i].sin_phi)*r;
        if (dop > Lx) {dop -= Lx;}
        if (dop < 0)  {dop += Lx;}
        xxx = dop*scale;
        dop = coordinates[i].y + dn2 * (-0.5*coordinates[i].sin_phi - 0.86602540378443864676372317075294*coordinates[i].cos_phi)*r;
        if (dop > Ly) {dop -= Ly;}
        if (dop < 0)  {dop += Ly;}
        yyy = dop*scale;
        xy_matrix[xxx][yyy]++;
    }
}
