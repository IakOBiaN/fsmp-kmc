class coord {
public:
   double x;
   double y;
   double z;
};
class ver {
public:
    double vol;
    coord place;
};
double interpol(ver v_000, ver v_010, ver v_100, ver v_110, ver v_001, ver v_011, ver v_101, ver v_111, coord v_x) {
                double d_x, d_y, d_z, vol_inter;
                d_x=(v_x.x-v_001.place.x)/(v_101.place.x-v_001.place.x);
                d_y=(v_x.y-v_000.place.y)/(v_010.place.y-v_000.place.y);
                d_z=(v_x.z-v_000.place.z)/(v_001.place.z-v_000.place.z);
                vol_inter=((v_000.vol*(1-d_x)+v_100.vol*d_x)*(1-d_y)+(v_010.vol*(1-d_x)+v_110.vol*d_x)*d_y)*(1-d_z)+((v_001.vol*(1-d_x)+v_101.vol*d_x)*(1-d_y)+(v_011.vol*(1-d_x)+v_111.vol*d_x)*d_y)*d_z;
                return vol_inter;
                }
