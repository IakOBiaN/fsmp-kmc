// Effective molecular area computed from the numerical potential, in A^2.
//
// For every approach direction the contact distance r0 is found as the inner
// boundary of the repulsive core: going outward from min_dist, the first
// distance where the pair energy stops being positive (linearly interpolated
// between grid points), minimized over all orientations of the second
// molecule. r0 is a center-to-center distance, i.e. an effective diameter, so
// the molecule itself contributes r0/2 in each direction; the returned value
// is the area of that halved contour,
//     A = 1/2 * Integral (r0(alpha)/2)^2 d(alpha),
// which for a spherical core (r0 = min_dist everywhere) reduces exactly to the
// old min_dist formula  pi * min_dist^2 / 4.
double molecule_area()
{
	int n_per = (int)llround(ff_fold_deg / da);   // stored period without the duplicate endpoint
	int n_r = cut_index + 1;                      // scan up to max_dist
	double r_far = min_dist + (n_r - 1) * dr;
	double sum_r0_2 = 0;
	for (int j1 = 0; j1 < n_per; j1++)
	{
		double r0_min = r_far;
		for (int j2 = 0; j2 < ff_nang; j2++)
		{
			double prev = FF(0, j1, j2);
			double r0 = min_dist;                   // no positive wall at min_dist: the core ends here
			if (prev > 0)
			{
				r0 = r_far;                           // stays there if the energy never crosses zero
				for (int i = 1; i < n_r; i++)
				{
					double cur = FF(i, j1, j2);
					if (cur <= 0)
					{
						r0 = min_dist + (i - 1) * dr + dr * prev / (prev - cur);
						break;
					}
					prev = cur;
				}
			}
			if (r0 < r0_min) { r0_min = r0; }
		}
		sum_r0_2 += r0_min * r0_min;
	}
	double copies = 360.0 / ff_fold_deg;          // a folded grid stores one period of the sweep
	return 0.125 * copies * sum_r0_2 * (da * PI / 180.0);
}
