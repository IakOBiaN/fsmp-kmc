#!/usr/bin/env python3
"""Compare the deterministic initial energy printed by the engine with a pin.

Reads the run log, takes the first "Energy:" value (printed for the central
cell before any Monte Carlo step) and checks it against the expected value.

Usage: check_energy.py <run.log> <expected_kJmol> <tolerance_kJmol>
"""
import re
import sys

log, expected, tol = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])

with open(log) as f:
    for line in f:
        m = re.search(r"Energy:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", line)
        if m:
            value = float(m.group(1))
            if abs(value - expected) <= tol:
                print("OK  initial energy %g kJ/mol (expected %g +- %g)"
                      % (value, expected, tol))
                sys.exit(0)
            print("FAIL initial energy %g kJ/mol, expected %g +- %g"
                  % (value, expected, tol))
            sys.exit(1)

print("FAIL: no 'Energy:' line found in %s" % log)
sys.exit(1)
