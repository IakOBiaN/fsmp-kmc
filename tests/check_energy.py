#!/usr/bin/env python3
"""Compare a deterministic energy printed by the engine with a pin.

Reads the run log and takes the value after the first occurrence of the label
("Energy:" by default, i.e. the central-cell energy printed before any Monte
Carlo step; pass another label to pin a different line, e.g. the optimizer's
"Final energy per molecule:").

Usage: check_energy.py <run.log> <expected_kJmol> <tolerance_kJmol> [label]
"""
import re
import sys

log, expected, tol = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
label = sys.argv[4] if len(sys.argv) > 4 else "Energy:"

with open(log) as f:
    for line in f:
        m = re.search(re.escape(label) + r"\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", line)
        if m:
            value = float(m.group(1))
            if abs(value - expected) <= tol:
                print("OK  initial energy %g kJ/mol (expected %g +- %g)"
                      % (value, expected, tol))
                sys.exit(0)
            print("FAIL initial energy %g kJ/mol, expected %g +- %g"
                  % (value, expected, tol))
            sys.exit(1)

print("FAIL: no '%s' line found in %s" % (label, log))
sys.exit(1)
