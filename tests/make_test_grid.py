#!/usr/bin/env python3
"""Make the small committed test grid by striding an existing v2 binary.

Keeps every N-th point along the distance and both angles, and scales dr/da in
the header accordingly. The angular grid must keep its endpoint, so N has to
divide n_ang-1. Assumes a little-endian host (the file format is little-endian).

The committed tests/data/TMA_simple_2020_s4.v2.bin was produced with:
    python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \
            tests/data/TMA_simple_2020_s4.v2.bin 4

Usage: make_test_grid.py <input.v2.bin> <output.v2.bin> <stride>
"""
import struct
import sys
from array import array

inp, outp, stride = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(inp, "rb") as f:
    hdr = bytearray(f.read(64))
    payload = f.read()

magic, ver, dtype, nd, na = struct.unpack_from("<4sIIII", hdr, 0)
min_d, dr, da, fold = struct.unpack_from("<dddd", hdr, 20)
assert magic == b"FSMP", "not an FSMP forcefield"
assert ver == 2, "format version %d, expected 2" % ver
assert stride >= 1, "stride must be >= 1"
assert (na - 1) % stride == 0, "stride %d does not divide n_ang-1 = %d" % (stride, na - 1)

fmt = "f" if dtype == 1 else "d"
vals = array(fmt)
vals.frombytes(payload)
assert len(vals) == nd * na * na, "payload size mismatch"

na2 = (na - 1) // stride + 1
nd2 = (nd - 1) // stride + 1
out = array(fmt)
for i in range(0, nd, stride):
    for j in range(0, na, stride):
        base = (i * na + j) * na
        out.extend(vals[base + k] for k in range(0, na, stride))

struct.pack_into("<II", hdr, 12, nd2, na2)
struct.pack_into("<dd", hdr, 28, dr * stride, da * stride)
with open(outp, "wb") as f:
    f.write(hdr)
    f.write(out.tobytes())

print("wrote %s: n_dist %d->%d  n_ang %d->%d  dr=%g da=%g fold=%g dtype=%s"
      % (outp, nd, nd2, na, na2, dr * stride, da * stride, fold,
         "float" if dtype == 1 else "double"))
