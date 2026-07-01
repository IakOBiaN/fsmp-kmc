#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>

// Reads a compact binary forcefield produced by tools/pack_forcefield. The grid is a
// flat block of energies (J/mol) in (distance, angle1, angle2) order, preceded by a
// 64-byte header. One extra distance-row of zeros is appended as a guard so the
// trilinear interpolation can read the (r + dr) neighbour at r = max_dist in range.
void read_forcefield(const char * filename, vector<double> & ff, int & n_ang,
                     double & min_dist, double & max_dist, double & dr, double & da,
                     double & fold_deg)
{
	FILE * pf = fopen(filename, "rb");
	if (!pf) { cerr << "ERROR: cannot open forcefield \"" << filename << "\"" << endl; exit(1); }

	char magic[4];
	uint32_t version = 0, dtype = 0, n_dist_u = 0, n_ang_u = 0;
	bool ok = true;
	ok &= fread(magic, 1, 4, pf) == 4;
	ok &= fread(&version,  4, 1, pf) == 1;
	ok &= fread(&dtype,    4, 1, pf) == 1;
	ok &= fread(&n_dist_u, 4, 1, pf) == 1;
	ok &= fread(&n_ang_u,  4, 1, pf) == 1;
	ok &= fread(&min_dist, 8, 1, pf) == 1;
	ok &= fread(&dr,       8, 1, pf) == 1;
	ok &= fread(&da,       8, 1, pf) == 1;
	ok &= fread(&fold_deg, 8, 1, pf) == 1;
	char skip[64];
	fread(skip, 1, 64 - 52, pf);   // header is padded to 64 bytes
	if (!ok || memcmp(magic, "FSMP", 4) != 0)
	{
		cerr << "ERROR: \"" << filename << "\" is not an FSMP forcefield. "
		     << "Convert the ASCII potential with tools/pack_forcefield first." << endl;
		exit(1);
	}
	if (version != 1) { cerr << "ERROR: unsupported forcefield version " << version << endl; exit(1); }

	int n_dist = (int)n_dist_u;
	n_ang = (int)n_ang_u;
	max_dist = min_dist + (n_dist - 1) * dr;

	size_t slab = (size_t)n_ang * n_ang;
	size_t count = (size_t)n_dist * slab;
	ff.assign((size_t)(n_dist + 1) * slab, 0.0);   // last distance-row stays zero (guard)

	bool read_ok = false;
	if (dtype == 0)
	{
		read_ok = fread(ff.data(), sizeof(double), count, pf) == count;
	}
	else if (dtype == 1)
	{
		vector<float> tmp(count);
		read_ok = fread(tmp.data(), sizeof(float), count, pf) == count;
		for (size_t m = 0; m < count; ++m) ff[m] = tmp[m];
	}
	else { cerr << "ERROR: unknown forcefield dtype " << dtype << endl; exit(1); }
	fclose(pf);
	if (!read_ok) { cerr << "ERROR: forcefield \"" << filename << "\" is truncated" << endl; exit(1); }

	cout << "forcefield: n_dist=" << n_dist << " n_ang=" << n_ang
	     << " r=[" << min_dist << ", " << max_dist << "] dr=" << dr
	     << " da=" << da << " fold=" << fold_deg << " deg" << endl;
}
