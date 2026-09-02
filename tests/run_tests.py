#!/usr/bin/env python3
"""Regression tests. Run from anywhere with any Python 3:

    python3 tests/run_tests.py          (Linux, macOS)
    python tests\\run_tests.py           (Windows, MinGW g++ on PATH)

1. pack_forcefield round-trip on a synthetic grid (no data needed).
2. The engine reports the version baked in from version.h (--version).
3. The full engine on the small committed grid in samples/potentials/.
4. The full engine on the real TMA simple potential; skipped when the
   potential is not present (it is a separate download, see README).
5. The unit-cell optimizer ("calculate") on the small committed grid,
   seeded for determinism; pins the converged energy.
6. The same optimization started from a cell with hard-core overlaps: the
   scaling stage must grow it out of the overlap and reach the same optimum.
7. The quickstart: the configuration a first-time user runs, exactly as
   shipped, on the demonstration potential that ships with it. Half a
   minute, and it pins the initial energy of its cell.

The engine is built once with -ffp-contract=off so the optimizer trajectory
(and its pin) is reproducible across compilers. Override the compiler with
the CXX environment variable (default g++).

samples/potentials/TMA_simple_2020_coarse_demo.v2.bin is a coarse (every 4th
point) copy of the distributed TMA simple binary, small enough to live in
git. It is what the quickstart runs on, so the suite exercises the same file
a first-time user does. Regenerable with:
    python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \\
            samples/potentials/TMA_simple_2020_coarse_demo.v2.bin 4
"""

import os
import subprocess
import sys
from pathlib import Path

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
BUILD = TESTS / "build"
EXE = ".exe" if os.name == "nt" else ""
CXX = os.environ.get("CXX", "g++")


def run(cmd, cwd=TESTS, **kw):
    """Run a command (from tests/ by default), aborting the suite when it
    fails."""
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, **kw)
    if result.returncode != 0:
        sys.exit(result.returncode)


def compile_cpp(source, out, *flags):
    run([CXX, "-O2", *flags, source, "-o", out])


def engine(config, log_name, cwd=TESTS):
    log = BUILD / log_name
    with open(log, "w") as sink:
        run([BUILD / ("fsmp" + EXE), config], cwd=cwd, stdout=sink,
            stderr=subprocess.STDOUT)
    return log


def rejected(edits, expected_message):
    """Run the engine on hcp_small.txt with some keys replaced, and require
    that it refuses the file with a message naming the problem."""
    lines = (TESTS / "hcp_small.txt").read_text().splitlines()
    for key, value in edits.items():
        for i, line in enumerate(lines):
            if line.split("=")[0].strip() == key:
                lines[i] = f"{key} = {value}"
                break
        else:
            sys.exit(f"test bug: hcp_small.txt has no key {key}")
    bad = BUILD / "invalid.txt"
    bad.write_text("\n".join(lines) + "\n")

    result = subprocess.run([str(BUILD / ("fsmp" + EXE)), "build/invalid.txt"],
                            cwd=TESTS, capture_output=True, text=True)
    output = result.stdout + result.stderr
    edited = ", ".join(f"{k} = {v}" for k, v in edits.items())
    if result.returncode == 0:
        sys.exit(f"FAIL: the engine accepted {edited}")
    if expected_message not in output:
        sys.exit(f"FAIL: {edited} was refused, but the message does not mention "
                 f"{expected_message!r}:\n{output.strip()}")
    print(f"OK  refused {edited}")


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

print("== [1/8] pack_forcefield round-trip on a synthetic grid ==", flush=True)
compile_cpp(TESTS.parent / "tools" / "pack_forcefield.cpp",
            BUILD / ("pack" + EXE), "-Wall", "-Wextra")
run([sys.executable, "test_pack_roundtrip.py", BUILD / ("pack" + EXE), BUILD])

compile_cpp(TESTS.parent / "fsmp.cpp", BUILD / ("fsmp" + EXE),
            "-ffp-contract=off")

print("== [2/8] the engine reports its version ==", flush=True)
version = subprocess.run([str(BUILD / ("fsmp" + EXE)), "--version"],
                         cwd=TESTS, capture_output=True, text=True)
if version.returncode != 0 or not version.stdout.startswith("FSMP-kMC "):
    sys.exit(f"--version failed: {version.stdout}{version.stderr}")
print(version.stdout.strip())

print("== [3/8] engine on the small committed grid ==", flush=True)
pin(engine("hcp_small.txt", "hcp_small.log"), -61.7914, 0.001)

print("== [4/8] engine on the full TMA simple potential ==", flush=True)
if (TESTS.parent / "forcefields" / "TMA_simple_2020.v2.bin").is_file():
    pin(engine("hcp_full.txt", "hcp_full.log"), -62.8777, 0.001)
else:
    print("SKIP: forcefields/TMA_simple_2020.v2.bin not present")

print("== [5/8] unit-cell optimizer on the small committed grid ==", flush=True)
pin(engine("optimize_small.txt", "optimize_small.log"), -62.211, 0.05,
    "Final energy per molecule:")

print("== [6/8] unit-cell optimizer from an overlapping start ==", flush=True)
pin(engine("optimize_overlap.txt", "optimize_overlap.log"), -62.211, 0.05,
    "Final energy per molecule:")

# The quickstart is the first thing a new user runs, so the suite runs it
# too: the shipped configuration, from the repository root as its own
# comment says, on the demonstration grid that ships next to it. Its output
# files land in the working directory, so they are cleared around the run.
print("== [7/8] the quickstart configuration, as shipped ==", flush=True)
QUICKSTART_OUTPUT = ("quickstart_0_unit_cell.xyz", "quickstart_1_trajectory.xyz",
                     "quickstart_2_statistics.dat")
for name in QUICKSTART_OUTPUT:
    (REPO / name).unlink(missing_ok=True)
pin(engine(Path("configs") / "tma_quickstart_demo.txt", "quickstart.log",
           cwd=REPO), -61.7304, 0.001)
for name in QUICKSTART_OUTPUT:
    (REPO / name).unlink(missing_ok=True)


# A parameter file that parses but cannot be run with must stop the program
# at once, naming the key and the line, instead of producing NaNs or looping
# forever hours into a run.
print("== [8/8] the engine refuses impossible parameters ==", flush=True)
rejected({"temp_from": "nan"}, "is not a finite number")
rejected({"um_to": "inf"}, "is not a finite number")
rejected({"nSteps": "3000000000"}, "does not fit in a 32-bit integer")
rejected({"temp_to": "400", "temp_step": "0"}, "would never reach the end")
rejected({"free_space": "0.5"}, "must be in [0, 0.5)")
rejected({"nStepsEq": "3"}, "must not exceed nSteps")
rejected({"uc_in_x": "0"}, "must be at least one unit cell")
rejected({"delta": "0"}, "must be a positive maximal displacement")
rejected({"temperature_in_transition_zone": "0"}, "must be a positive temperature")
# nSteps itself fits in an int, but nSteps x 264 particles does not
rejected({"nSteps": "10000000"}, "iterations, which does not fit")

print("ALL TESTS PASSED")
