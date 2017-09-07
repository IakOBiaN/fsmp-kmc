void layer_map(int &nPart, vector <state> &coordinates, vector <vector <double>> &xy_matrix)
{
    double dn2 = 0.5;

    for(int i = 0; i < nPart; i++)
    {
        double r = (1 - 2 * RanGen.Random());
        int xxx = (coordinates[i].x + dn2 * sin(coordinates[i].tetta)*cos(coordinates[i].phi)*r + 0.5)*10;
        int yyy = (coordinates[i].y + dn2 * sin(coordinates[i].tetta)*sin(coordinates[i].phi)*r + 0.5)*10;
        xy_matrix[xxx][yyy]++;
    }
}
