#!/bin/bash
# Regression tests. Run from anywhere:  ./tests/run_tests.sh
# Needs a C++ compiler (override with CXX=...) and python3.
#
# 1. pack_forcefield round-trip on a synthetic grid (no data needed).
# 2. The full engine on the small committed grid in tests/data/.
# 3. The full engine on the real TMA simple potential; skipped when the
#    potential is not present (it is a separate download, see README).
# 4. The unit-cell optimizer ("calculate") on the small committed grid,
#    seeded for determinism; pins the converged energy.
#
# The engine is built once with -ffp-contract=off so the optimizer trajectory
# (and its pin) is reproducible across compilers.
#
# tests/data/TMA_simple_2020_s4.v2.bin is a coarse (every 4th point) copy of
# the distributed TMA simple binary, regenerable with:
#   python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \
#           tests/data/TMA_simple_2020_s4.v2.bin 4
set -e
cd "$(dirname "$0")"
CXX="${CXX:-g++}"
mkdir -p build
# the engine never overwrites existing outputs, so clear the previous run's
# files to keep the output names canonical
rm -f 0_*.xyz 1_*.xyz 2_*.dat

echo "== [1/4] pack_forcefield round-trip on a synthetic grid =="
"$CXX" -O2 -Wall -Wextra ../tools/pack_forcefield.cpp -o build/pack
python3 test_pack_roundtrip.py build/pack build

"$CXX" -O2 -ffp-contract=off ../fsmp.cpp -o build/fsmp

echo "== [2/4] engine on the small committed grid =="
./build/fsmp hcp_small.txt > build/hcp_small.log
python3 check_energy.py build/hcp_small.log -61.7449 0.001

echo "== [3/4] engine on the full TMA simple potential =="
if [ -f ../forcefields/TMA_simple_2020.v2.bin ]; then
    ./build/fsmp hcp_full.txt > build/hcp_full.log
    python3 check_energy.py build/hcp_full.log -62.8605 0.001
else
    echo "SKIP: forcefields/TMA_simple_2020.v2.bin not present"
fi

echo "== [4/4] unit-cell optimizer on the small committed grid =="
./build/fsmp optimize_small.txt > build/optimize_small.log
python3 check_energy.py build/optimize_small.log -62.2278 0.05 "Final energy per molecule:"

echo "ALL TESTS PASSED"
