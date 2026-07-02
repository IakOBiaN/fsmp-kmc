#!/bin/bash
# Regression tests. Run from anywhere:  ./tests/run_tests.sh
# Needs a C++ compiler (override with CXX=...) and python3.
#
# 1. pack_forcefield round-trip on a synthetic grid (no data needed).
# 2. The full engine on the small committed grid in tests/data/.
# 3. The full engine on the real TMA simple potential; skipped when the
#    potential is not present (it is a separate download, see README).
#
# tests/data/TMA_simple_2020_s4.v2.bin is a coarse (every 4th point) copy of
# the distributed TMA simple binary, regenerable with:
#   python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \
#           tests/data/TMA_simple_2020_s4.v2.bin 4
set -e
cd "$(dirname "$0")"
CXX="${CXX:-g++}"
mkdir -p build

echo "== [1/3] pack_forcefield round-trip on a synthetic grid =="
"$CXX" -O2 -Wall -Wextra ../tools/pack_forcefield.cpp -o build/pack
python3 test_pack_roundtrip.py build/pack build

echo "== [2/3] engine on the small committed grid =="
"$CXX" -O3 hcp_small.cpp -o build/hcp_small
./build/hcp_small > build/hcp_small.log
python3 check_energy.py build/hcp_small.log -61.7449 0.001

echo "== [3/3] engine on the full TMA simple potential =="
if [ -f ../forcefields/TMA_simple_2020.v2.bin ]; then
    "$CXX" -O3 hcp_full.cpp -o build/hcp_full
    ./build/hcp_full > build/hcp_full.log
    python3 check_energy.py build/hcp_full.log -62.8605 0.001
else
    echo "SKIP: forcefields/TMA_simple_2020.v2.bin not present"
fi

echo "ALL TESTS PASSED"
