#!/usr/bin/env python3
"""Round-trip test of tools/pack_forcefield.

Generates a small synthetic ASCII grid (120-deg periodic, exchange-symmetric),
packs it in several modes, and verifies every byte of the output against an
independent re-implementation of the packing rules: kcal/mol -> J/mol scaling,
molecule-exchange symmetrization, rotational folding, striding, float rounding.
Also checks that the tool refuses invalid inputs (bad fold period, truncated
file) and cleans up after itself.

Usage: test_pack_roundtrip.py <pack-binary> <workdir>
"""
import math
import os
import struct
import subprocess
import sys

PACK = os.path.abspath(sys.argv[1])
WORK = os.path.abspath(sys.argv[2])

MIN_DIST, DR, N_DIST = 5.0, 0.5, 4
DA, N_ANG = 30.0, 13  # angles 0..360 inclusive


def energy(r, t1, t2):
    return (1.0 / r) * math.cos(3 * math.radians(t1)) * math.cos(3 * math.radians(t2)) + 0.01 * r


def write_ascii(path):
    """Write the grid and return the values exactly as the C tool will parse them."""
    grid = []  # grid[i][j][k], J/mol
    with open(path, "w") as f:
        for i in range(N_DIST):
            r = MIN_DIST + i * DR
            slab = []
            for j in range(N_ANG):
                row = []
                for k in range(N_ANG):
                    text = "%.12g" % energy(r, j * DA, k * DA)
                    f.write("%.6g %.6g %.6g %s\n" % (r, j * DA, k * DA, text))
                    row.append(float(text) * 4184.0)
                slab.append(row)
            grid.append(slab)
    return grid


def expected_payload(grid, fold_deg, stride):
    """Re-implement subsample -> symmetrize -> fold on the parsed grid."""
    sub = [[[grid[i][j * stride][k * stride]
             for k in range((N_ANG - 1) // stride + 1)]
            for j in range((N_ANG - 1) // stride + 1)]
           for i in range(0, N_DIST, stride)]
    n_ang = len(sub[0])
    da = DA * stride
    half, full = round(180.0 / da), round(360.0 / da)
    for s in sub:  # exchange symmetry, same in-place order as the C tool
        for j in range(n_ang):
            jp = j + half
            if jp >= full:
                jp -= full
            for k in range(n_ang):
                kp = k + half
                if kp >= full:
                    kp -= full
                m = max(s[j][k], s[kp][jp])
                s[j][k] = m
                s[kp][jp] = m
    out_n_ang = n_ang
    if fold_deg < 359.999:
        out_n_ang = round(fold_deg / da) + 1
    flat = [s[j][k] for s in sub for j in range(out_n_ang) for k in range(out_n_ang)]
    header = (len(sub), out_n_ang, MIN_DIST, DR * stride, da,
              fold_deg if fold_deg < 359.999 else 360.0)
    return header, flat


def read_bin(path):
    with open(path, "rb") as f:
        hdr = f.read(64)
        payload = f.read()
    magic, ver, dtype, nd, na = struct.unpack_from("<4sIIII", hdr, 0)
    min_d, dr, da, fold = struct.unpack_from("<dddd", hdr, 20)
    assert magic == b"FSMP", "bad magic %r" % magic
    assert ver == 2, "format version %d, expected 2" % ver
    count = len(payload) // (4 if dtype == 1 else 8)
    vals = struct.unpack("<%d%s" % (count, "f" if dtype == 1 else "d"), payload)
    return dtype, (nd, na, min_d, dr, da, fold), list(vals)


def run_pack(args):
    return subprocess.run([PACK] + args, cwd=WORK,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


def check_case(grid, name, use_float, fold_deg, stride):
    out = "toy_%s.bin" % name
    args = (["--float"] if use_float else []) \
        + (["--stride", str(stride)] if stride > 1 else []) \
        + ["toy.dat", out] + ([str(fold_deg)] if fold_deg < 360 else [])
    rc = run_pack(args)
    assert rc == 0, "%s: pack failed (rc=%d)" % (name, rc)
    dtype, header, vals = read_bin(os.path.join(WORK, out))
    assert dtype == (1 if use_float else 0), "%s: dtype %d" % (name, dtype)
    exp_header, exp = expected_payload(grid, fold_deg, stride)
    got_geom = (header[0], header[1])
    assert (exp_header[0], exp_header[1]) == got_geom, \
        "%s: geometry %s vs expected %s" % (name, got_geom, exp_header[:2])
    for h, e, what in zip(header[2:], exp_header[2:], ("min_dist", "dr", "da", "fold")):
        assert abs(h - e) < 1e-12, "%s: header %s = %g, expected %g" % (name, what, h, e)
    if use_float:
        exp = [struct.unpack("<f", struct.pack("<f", v))[0] for v in exp]
    assert len(vals) == len(exp), "%s: %d values, expected %d" % (name, len(vals), len(exp))
    bad = sum(1 for a, b in zip(vals, exp) if a != b)
    assert bad == 0, "%s: %d/%d payload values differ" % (name, bad, len(exp))
    print("OK  %-14s dtype=%-6s n_dist=%d n_ang=%d fold=%g stride=%d"
          % (name, "float" if use_float else "double",
             header[0], header[1], header[5], stride))


def main():
    os.makedirs(WORK, exist_ok=True)
    grid = write_ascii(os.path.join(WORK, "toy.dat"))

    check_case(grid, "d",       False, 360, 1)
    check_case(grid, "f",       True,  360, 1)
    check_case(grid, "d_fold",  False, 120, 1)
    check_case(grid, "f_fold",  True,  120, 1)
    check_case(grid, "d_s2",    False, 360, 2)
    check_case(grid, "f_fold_s2", True, 120, 2)

    # the toy grid is 120-deg periodic but NOT 180-deg periodic: folding at 180
    # must be refused with the default tolerance, and no output may be left over
    rc = run_pack(["toy.dat", "toy_bad.bin", "180"])
    assert rc != 0, "fold at 180 deg unexpectedly succeeded"
    assert not os.path.exists(os.path.join(WORK, "toy_bad.bin")), "toy_bad.bin left behind"
    print("OK  refused a fold period the data does not have")

    # a truncated input must be rejected, not silently packed
    with open(os.path.join(WORK, "toy.dat")) as f:
        lines = f.readlines()
    with open(os.path.join(WORK, "toy_trunc.dat"), "w") as f:
        f.writelines(lines[: len(lines) - 5])
    rc = run_pack(["toy_trunc.dat", "toy_trunc.bin"])
    assert rc != 0, "truncated input unexpectedly packed"
    assert not os.path.exists(os.path.join(WORK, "toy_trunc.bin")), "toy_trunc.bin left behind"
    print("OK  refused a truncated input")

    print("pack_forcefield round-trip: all cases passed")


if __name__ == "__main__":
    main()
