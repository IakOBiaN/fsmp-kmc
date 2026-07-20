# FSMP-kMC build helper.
#
#   make            build the simulation program (fsmp.out)
#   make pack       build the forcefield converter (pack.out)
#   make windows    build native Windows binaries (fsmp.exe, pack.exe);
#                   needs a MinGW g++ (w64devkit or MSYS2) on PATH
#   make test       run the regression test suite (needs a Python 3; on
#                   Windows pass PYTHON=python or a full interpreter path)
#   make bundle     assemble a release-layout bundle in dist/ (Windows;
#                   elsewhere run tools/make_bundle.py --build-engine)
#   make clean      remove built binaries
#
# Run the program with a parameter file, from the repository root:
#   ./fsmp.out configs/tma_acid_hcp.txt
# Paths in the parameter file are relative to the working directory and all
# output files are written there.

CXX      ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra
PYTHON   ?= python3

.PHONY: all pack windows test bundle clean

all: fsmp.out

fsmp.out: fsmp.cpp program_body.cpp includes.h $(wildcard *.h)
	$(CXX) $(CXXFLAGS) fsmp.cpp -o $@

pack: pack.out

pack.out: tools/pack_forcefield.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

# -static bakes the MinGW runtime in, so the binaries run on any Windows
# machine with no extra DLLs (the release workflow builds through this
# same target)
windows: fsmp.exe pack.exe

# Version resources (icon and the Explorer Properties fields) for the
# native binaries. The numeric a,b,c,0 quad is derived here from the
# "a.b.c" string in version.h, the single project-version source.
comma := ,
RC_VERSION = $(subst .,$(comma),$(shell sed -n 's/.*FSMP_VERSION "\([^"]*\)".*/\1/p' version.h)),0

fsmp.res.o: tools/fsmp.rc version.h gui/studio.ico
	windres -DFSMP_RC_VERSION=$(RC_VERSION) --include-dir . $< $@

pack.res.o: tools/pack.rc version.h gui/studio.ico
	windres -DFSMP_RC_VERSION=$(RC_VERSION) --include-dir . $< $@

fsmp.exe: fsmp.cpp program_body.cpp includes.h $(wildcard *.h) fsmp.res.o
	$(CXX) $(CXXFLAGS) -static fsmp.cpp fsmp.res.o -o $@

pack.exe: tools/pack_forcefield.cpp pack.res.o
	$(CXX) $(CXXFLAGS) -static $< pack.res.o -o $@

test:
	$(PYTHON) tests/run_tests.py

bundle: windows
	$(PYTHON) tools/make_bundle.py

clean:
	rm -f fsmp.out pack.out *.res.o configs/*.out
	rm -rf tests/build
