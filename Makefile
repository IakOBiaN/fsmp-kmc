# FSMP-kMC build helper.
#
#   make            build the simulation program (fsmp.out)
#   make pack       build the forcefield converter (pack.out)
#   make test       run the regression test suite
#   make clean      remove built binaries
#
# Run the program with a parameter file, from the repository root:
#   ./fsmp.out configs/tma_acid_hcp.txt
# Paths in the parameter file are relative to the working directory and all
# output files are written there.

CXX      ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra

.PHONY: all pack test clean

all: fsmp.out

fsmp.out: fsmp.cpp program_body.cpp includes.h $(wildcard *.h)
	$(CXX) $(CXXFLAGS) fsmp.cpp -o $@

pack: pack.out

pack.out: tools/pack_forcefield.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test:
	bash tests/run_tests.sh

clean:
	rm -f fsmp.out pack.out configs/*.out
	rm -rf tests/build
