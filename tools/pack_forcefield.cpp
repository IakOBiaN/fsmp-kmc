// pack_forcefield: convert a legacy ASCII forcefield into the compact binary format.
//
// Legacy files store one grid point per line as four ASCII numbers
//     r(A)   angle1(deg)   angle2(deg)   energy(kcal/mol)
// on a regular grid (r slowest, angle1 next, angle2 fastest). The coordinates are
// redundant, so the binary keeps only the energy values (J/mol) in grid order behind a
// 64-byte header that records the geometry. The run time reads this format; ASCII is no
// longer read at run time. Use this tool once to convert a downloaded or self-generated
// potential.
//
// Bit-faithful to the old loader: energies are scaled kcal/mol -> J/mol (x4184), and the
// molecule-exchange symmetry U(r,t1,t2) = U(r,t2+180,t1+180) is enforced (max of the two
// partners).
//
// Optional rotational folding: pass a period in degrees (e.g. 120 for the C3 trimesic
// acid, 180 for a C2 molecule) to store a single angular period. Periodicity is checked
// only in the physical region (where the energy is below the runtime hard-core cap at
// T_ref; above it all values fold losslessly to the cap), and the tool refuses to fold if
// the physical deviation exceeds the tolerance. With no period (or 360) the grid is stored
// in full, and the tool reports what a 120 / 180 deg fold would cost.
//
// Streams one distance slab (n_ang x n_ang) at a time, so it converts arbitrarily large
// grids in O(n_ang^2) memory.
//
// Usage:  pack_forcefield <input.dat> <output.bin> [fold_deg] [tol_Jmol] [T_ref_K]
// Build:  clang++ -O3 tools/pack_forcefield.cpp -o pack
//
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <vector>

static const char     MAGIC[4]       = {'F', 'S', 'M', 'P'};
static const uint32_t FORMAT_VERSION = 1;
static const int      HEADER_BYTES   = 64;
static const double   E_INF          = 75.0;         // hard-core cap, in units of kT
static const double   R              = 8.314462618;  // gas constant, J/(mol K)

// Largest deviation between points one period (p steps) apart, along both angular axes of
// a slab, counting only pairs where at least one partner is below thr (the physical
// region: pairs where both are above the runtime cap fold losslessly). -1 if p is invalid.
static double slab_period_dev(const std::vector<double>& s, int n_ang, int p, double thr) {
    if (p <= 0 || p >= n_ang) return -1.0;
    double maxd = 0.0;
    for (int j = 0; j < n_ang; ++j)
        for (int k = 0; k < n_ang; ++k) {
            double v = s[(size_t)j * n_ang + k];
            if (j + p < n_ang) { double w = s[(size_t)(j + p) * n_ang + k]; if ((v < thr || w < thr) && fabs(v - w) > maxd) maxd = fabs(v - w); }
            if (k + p < n_ang) { double w = s[(size_t)j * n_ang + (k + p)]; if ((v < thr || w < thr) && fabs(v - w) > maxd) maxd = fabs(v - w); }
        }
    return maxd;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <input.dat> <output.bin> [fold_deg] [tol_Jmol] [T_ref_K]\n", argv[0]); return 2; }
    const char* in_path  = argv[1];
    const char* out_path = argv[2];
    double      fold_deg = (argc > 3) ? atof(argv[3]) : 360.0;
    double      fold_tol = (argc > 4) ? atof(argv[4]) : 100.0;   // max physical |dU|, J/mol
    double      T_ref    = (argc > 5) ? atof(argv[5]) : 1000.0;  // reference temperature, K
    double      thr      = E_INF * R * T_ref;                    // runtime cap at T_ref, J/mol

    FILE* fin = fopen(in_path, "r");
    if (!fin) { fprintf(stderr, "cannot open %s\n", in_path); return 1; }

    // --- detect geometry from the first distance block (reads one block, not the file) ---
    double first_dist = 0, prev_a1 = 0, da = 0, dr = 0, r, a1, a2, e;
    bool have_da = false; long block_lines = 1;
    if (fscanf(fin, "%lf %lf %lf %lf", &r, &a1, &a2, &e) != 4) { fprintf(stderr, "empty input\n"); return 1; }
    first_dist = r; prev_a1 = a1;
    while (fscanf(fin, "%lf %lf %lf %lf", &r, &a1, &a2, &e) == 4) {
        if (r != first_dist) { dr = r - first_dist; break; }
        if (!have_da && a1 != prev_a1) { da = a1 - prev_a1; have_da = true; }
        prev_a1 = a1; ++block_lines;
    }
    long n_ang_l = (long)llround(sqrt((double)block_lines));
    if (n_ang_l * n_ang_l != block_lines) { fprintf(stderr, "first block (%ld lines) is not a perfect square\n", block_lines); return 1; }
    if (!have_da || dr == 0) { fprintf(stderr, "could not detect da/dr\n"); return 1; }
    int n_ang = (int)n_ang_l;
    int half = (int)llround(180.0 / da), full = (int)llround(360.0 / da);

    int    fold_p = 0, out_n_ang = n_ang; double out_fold = 360.0;
    bool   do_fold = (fold_deg < 359.999);
    if (do_fold) {
        fold_p = (int)llround(fold_deg / da);
        if (fold_p <= 0 || fabs(fold_p * da - fold_deg) > 1e-9) { fprintf(stderr, "fold_deg not a multiple of da=%.4f\n", da); return 1; }
        out_n_ang = fold_p + 1; out_fold = fold_deg;
    }

    // --- write a placeholder header (n_dist is patched once we have counted slabs) ---
    FILE* fout = fopen(out_path, "wb");
    if (!fout) { fprintf(stderr, "cannot open %s\n", out_path); return 1; }
    uint32_t dtype = 0, nd0 = 0, na = (uint32_t)out_n_ang;
    long hb = 0;
    hb += fwrite(MAGIC, 1, 4, fout);
    hb += fwrite(&FORMAT_VERSION, 4, 1, fout) * 4;
    hb += fwrite(&dtype, 4, 1, fout) * 4;
    long off_n_dist = hb;
    hb += fwrite(&nd0, 4, 1, fout) * 4;
    hb += fwrite(&na, 4, 1, fout) * 4;
    hb += fwrite(&first_dist, 8, 1, fout) * 8;
    hb += fwrite(&dr, 8, 1, fout) * 8;
    hb += fwrite(&da, 8, 1, fout) * 8;
    hb += fwrite(&out_fold, 8, 1, fout) * 8;
    char pad[HEADER_BYTES] = {0};
    fwrite(pad, 1, HEADER_BYTES - hb, fout);

    // --- stream every slab from the start of the file ---
    rewind(fin);
    size_t sn = (size_t)n_ang * n_ang;
    std::vector<double> slab(sn), out;
    if (do_fold) out.resize((size_t)out_n_ang * out_n_ang);
    int p120 = (int)llround(120.0 / da), p180 = (int)llround(180.0 / da);
    double dev120 = 0, dev180 = 0, devfold = 0;
    long n_dist = 0;
    while (true) {
        size_t got = 0;
        for (; got < sn; ++got) {
            if (fscanf(fin, "%lf %lf %lf %lf", &r, &a1, &a2, &e) != 4) break;
            slab[got] = e * 4184.0;
        }
        if (got == 0) break;                            // clean end of file
        if (got != sn) { fprintf(stderr, "truncated slab %ld (%zu/%zu lines)\n", n_dist, got, sn); fclose(fout); remove(out_path); return 1; }

        for (int j = 0; j < n_ang; ++j) { int jp = j + half; if (jp >= full) jp -= full;
            for (int k = 0; k < n_ang; ++k) { int kp = k + half; if (kp >= full) kp -= full;
                double& A = slab[(size_t)j * n_ang + k]; double& B = slab[(size_t)kp * n_ang + jp];
                double m = (A > B) ? A : B; A = B = m; } }

        double d1 = slab_period_dev(slab, n_ang, p120, thr); if (d1 > dev120) dev120 = d1;
        double d2 = slab_period_dev(slab, n_ang, p180, thr); if (d2 > dev180) dev180 = d2;

        const std::vector<double>* w = &slab;
        if (do_fold) {
            double df = slab_period_dev(slab, n_ang, fold_p, thr); if (df > devfold) devfold = df;
            if (df > fold_tol) { fprintf(stderr, "slab %ld: physical periodicity error %.6g J/mol > tol %.6g at %.0f deg; refusing to fold (raise the tolerance to force)\n", n_dist, df, fold_tol, fold_deg); fclose(fout); remove(out_path); return 1; }
            for (int j = 0; j < out_n_ang; ++j)
                for (int k = 0; k < out_n_ang; ++k)
                    out[(size_t)j * out_n_ang + k] = slab[(size_t)j * n_ang + k];
            w = &out;
        }
        fwrite(w->data(), sizeof(double), w->size(), fout);
        ++n_dist;
    }
    fclose(fin);
    if (n_dist == 0) { fprintf(stderr, "no slabs written\n"); fclose(fout); remove(out_path); return 1; }

    fseek(fout, off_n_dist, SEEK_SET);
    uint32_t nd = (uint32_t)n_dist;
    fwrite(&nd, 4, 1, fout);
    fclose(fout);

    double max_dist = first_dist + (n_dist - 1) * dr;
    printf("packed %s\n", out_path);
    printf("  n_dist=%ld n_ang=%d -> out_n_ang=%d  r=[%.4f, %.4f] dr=%.5f da=%.4f fold=%.0f\n",
           n_dist, n_ang, out_n_ang, first_dist, max_dist, dr, da, out_fold);
    printf("  physical periodicity (cap at %.0f K = %.4g J/mol): 120 deg = %.4g J/mol, 180 deg = %.4g J/mol\n",
           T_ref, thr, dev120, dev180);
    if (do_fold) printf("  folded at %.0f deg, physical |dU| = %.4g J/mol (tol %.4g)\n", fold_deg, devfold, fold_tol);
    double mb = (double)(HEADER_BYTES + (size_t)n_dist * out_n_ang * out_n_ang * 8) / (1024.0 * 1024.0);
    printf("  size: %.1f MB (dtype=double)\n", mb);
    return 0;
}
