void density_change(double &Lx, double &Ly, int &nPart, vector <state> &coordinates)
{
    double area_new,Lx_new,Ly_new;    // Linear sizes of the box after correction
		double ratio_x_to_y = Lx/Ly;

    area_new = nPart*1000.0/density/6.02214129;
		Lx_new = sqrt(area_new*ratio_x_to_y);
		Ly_new = Lx_new/ratio_x_to_y;

       for(int i = 0; i < nPart; i++)
          {
            // Coordinates of the molecules change proportionally to
            // the change of the box size
            coordinates[i].x = (Lx_new/Lx)*coordinates[i].x;
            coordinates[i].y = (Ly_new/Ly)*coordinates[i].y;
          }
       Lx = Lx_new;
       Ly = Ly_new;
}
