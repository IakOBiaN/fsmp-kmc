# FSMP-kMC build helper. The engine is built per configuration: a file in
# configs/ includes program_body.cpp and defines one simulation.
#
#   make CONFIG=configs/tma_acid_hcp.cpp   build a run binary next to the config
#   make pack                              build the forcefield converter (pack.out)
#   make test                              run the regression test suite
#   make clean                             remove built binaries
#
# The run binary is placed next to its config and must be started from the
# configs/ directory (potential paths in the configs are relative to it):
#   cd configs && ./tma_acid_hcp.out

CXX      ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra
CONFIG   ?= configs/tma_acid_hcp.cpp

BIN := $(CONFIG:.cpp=.out)

.PHONY: all pack test clean

all: $(BIN)

$(BIN): $(CONFIG) program_body.cpp includes.h $(wildcard *.h)
	$(CXX) $(CXXFLAGS) $< -o $@

pack: pack.out

pack.out: tools/pack_forcefield.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test:
	bash tests/run_tests.sh

clean:
	rm -f configs/*.out pack.out
	rm -rf tests/build
