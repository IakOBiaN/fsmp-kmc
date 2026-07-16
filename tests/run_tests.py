#!/usr/bin/env python3
"""Regression tests. Run from anywhere with any Python 3:

    python3 tests/run_tests.py          (Linux, macOS)
    python tests\\run_tests.py           (Windows, MinGW g++ on PATH)

1. pack_forcefield round-trip on a synthetic grid (no data needed).
2. The engine reports the version baked in from version.h (--version).
3. The full engine on the small committed grid in tests/data/.
4. The full engine on the real TMA simple potential; skipped when the
   potential is not present (it is a separate download, see README).
5. The unit-cell optimizer ("calculate") on the small committed grid,
   seeded for determinism; pins the converged energy.
6. The same optimization started from a cell with hard-core overlaps: the
   scaling stage must grow it out of the overlap and reach the same optimum.

The engine is built once with -ffp-contract=off so the optimizer trajectory
(and its pin) is reproducible across compilers. Override the compiler with
the CXX environment variable (default g++).

tests/data/TMA_simple_2020_s4.v2.bin is a coarse (every 4th point) copy of
the distributed TMA simple binary, regenerable with:
    python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \\
            tests/data/TMA_simple_2020_s4.v2.bin 4
"""

import os
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent
BUILD = TESTS / "build"
EXE = ".exe" if os.name == "nt" else ""
CXX = os.environ.get("CXX", "g++")


def run(cmd, **kw):
    """Run a command from tests/, aborting the suite when it fails."""
    result = subprocess.run([str(c) for c in cmd], cwd=TESTS, **kw)
    if result.returncode != 0:
        sys.exit(result.returncode)


def compile_cpp(source, out, *flags):
    run([CXX, "-O2", *flags, source, "-o", out])


def engine(config, log_name):
    log = BUILD / log_name
    with open(log, "w") as sink:
        run([BUILD / ("fsmp" + EXE), config], stdout=sink,
            stderr=subprocess.STDOUT)
    return log


def pin(log, expected, tolerance, label=None):
    cmd = [sys.executable, "check_energy.py", log, expected, tolerance]
    if label:
        cmd.append(label)
    run(cmd)


BUILD.mkdir(exist_ok=True)
# the engine never overwrites existing outputs, so clear the previous run's
# files to keep the output names canonical
for pattern in ("0_*.xyz", "1_*.xyz", "2_*.dat"):
    for stray in TESTS.glob(pattern):
        stray.unlink()

print("== [1/6] pack_forcefield round-trip on a synthetic grid ==", flush=True)
compile_cpp(TESTS.parent / "tools" / "pack_forcefield.cpp",
            BUILD / ("pack" + EXE), "-Wall", "-Wextra")
run([sys.executable, "test_pack_roundtrip.py", BUILD / ("pack" + EXE), BUILD])

compile_cpp(TESTS.parent / "fsmp.cpp", BUILD / ("fsmp" + EXE),
            "-ffp-contract=off")

print("== [2/6] the engine reports its version ==", flush=True)
version = subprocess.run([str(BUILD / ("fsmp" + EXE)), "--version"],
                         cwd=TESTS, capture_output=True, text=True)
if version.returncode != 0 or not version.stdout.startswith("FSMP-kMC "):
    sys.exit(f"--version failed: {version.stdout}{version.stderr}")
print(version.stdout.strip())

print("== [3/6] engine on the small committed grid ==", flush=True)
pin(engine("hcp_small.txt", "hcp_small.log"), -61.7449, 0.001)

print("== [4/6] engine on the full TMA simple potential ==", flush=True)
if (TESTS.parent / "forcefields" / "TMA_simple_2020.v2.bin").is_file():
    pin(engine("hcp_full.txt", "hcp_full.log"), -62.8605, 0.001)
else:
    print("SKIP: forcefields/TMA_simple_2020.v2.bin not present")

print("== [5/6] unit-cell optimizer on the small committed grid ==", flush=True)
pin(engine("optimize_small.txt", "optimize_small.log"), -62.2276, 0.05,
    "Final energy per molecule:")

print("== [6/6] unit-cell optimizer from an overlapping start ==", flush=True)
pin(engine("optimize_overlap.txt", "optimize_overlap.log"), -62.2276, 0.05,
    "Final energy per molecule:")

print("ALL TESTS PASSED")
