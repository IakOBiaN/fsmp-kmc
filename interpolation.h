class coord {
public:
   double x;
   double y;
   double z;
};
class ver {
public:
    double value;
    coord place;
};
double interpolation(double dist_n, double ang1, double ang2, double beta) {
                double d_x, d_y, d_z, energy_interpol;
                ver v_000, v_010, v_100, v_110, v_001, v_011, v_101, v_111;
                int temp_dist,temp_ang1,temp_ang2;
                coord v_find;
                temp_dist = (int)(dist_n);
                temp_ang1 = (int)(ang1 / da);
                temp_ang2 = (int)(ang2 / da);
                v_000.place.x = temp_dist;
                v_000.place.y = temp_ang1;
                v_000.place.z = temp_ang2;
                v_010.place.x = temp_dist;
                v_010.place.y = temp_ang1 + 1;
                v_010.place.z = temp_ang2;
                v_100.place.x = temp_dist + 1;
                v_100.place.y = temp_ang1;
                v_100.place.z = temp_ang2;
                v_110.place.x = temp_dist + 1;
                v_110.place.y = temp_ang1 + 1;
                v_110.place.z = temp_ang2;
                v_001.place.x = temp_dist;
                v_001.place.y = temp_ang1;
                v_001.place.z = temp_ang2 + 1;
                v_011.place.x = temp_dist;
                v_011.place.y = temp_ang1 + 1;
                v_011.place.z = temp_ang2 + 1;
                v_101.place.x = temp_dist + 1;
                v_101.place.y = temp_ang1;
                v_101.place.z = temp_ang2 + 1;
                v_111.place.x = temp_dist + 1;
                v_111.place.y = temp_ang1 + 1;
                v_111.place.z = temp_ang2 + 1;
                v_find.x = dist_n;
                v_find.y = ang1 / da;
                v_find.z = ang2 / da;
                v_000.value = forcefield[v_000.place.x][v_000.place.y][v_000.place.z];
                v_010.value = forcefield[v_010.place.x][v_010.place.y][v_010.place.z];
                v_100.value = forcefield[v_100.place.x][v_100.place.y][v_100.place.z];
                v_110.value = forcefield[v_110.place.x][v_110.place.y][v_110.place.z];
                v_001.value = forcefield[v_001.place.x][v_001.place.y][v_001.place.z];
                v_011.value = forcefield[v_011.place.x][v_011.place.y][v_011.place.z];
                v_101.value = forcefield[v_101.place.x][v_101.place.y][v_101.place.z];
                v_111.value = forcefield[v_111.place.x][v_111.place.y][v_111.place.z];
                d_x = (v_find.x - v_001.place.x) / (v_101.place.x - v_001.place.x);
                d_y = (v_find.y - v_000.place.y) / (v_010.place.y - v_000.place.y);
                d_z = (v_find.z - v_000.place.z) / (v_001.place.z - v_000.place.z);
                energy_interpol = ((v_000.value * (1 - d_x) + v_100.value * d_x) * (1 - d_y) + (v_010.value * (1 - d_x) + v_110.value * d_x) * d_y) * (1-d_z) + ((v_001.value * (1 - d_x) + v_101.value * d_x) * (1 - d_y) + (v_011.value * (1 - d_x) + v_111.value * d_x) * d_y) * d_z;
                if ((v_000.value > 0) && (v_010.value > 0) && (v_100.value > 0) && (v_110.value > 0) && (v_001.value > 0) && (v_011.value > 0) && (v_101.value > 0) && (v_111.value > 0) && (energy_interpol < 0))
                {
                  energy_interpol = -energy_interpol;
                }
                if (energy_interpol > (E_INF / beta))
                {
                  energy_interpol = E_INF / beta;
                }
                return energy_interpol;
                }
