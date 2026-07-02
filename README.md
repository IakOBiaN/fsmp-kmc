# FSMP-kMC

**Field-Stabilized Multiphase kinetic Monte Carlo.** A kinetic Monte Carlo
engine for the atomistic thermodynamics of rigid molecular crystals and
two-dimensional adsorption layers.

The code simulates two coexisting phases (a crystal and an ideal-gas reservoir)
in a single elongated cell. Two imposed inhomogeneous fields, a *damping* field
and an *external* field, stabilize this coexistence over a wide range of
temperature and pressure. This makes it possible to determine the free energy,
entropy and chemical potential of dense molecular layers from the equality of
chemical potentials in the coexisting phases.

## Method

This code accompanies the following study:

> S. S. Akimenko, V. A. Gorbunov and E. A. Ustinov, *Equilibrium structure of a
> dense trimesic acid monolayer on a homogeneous solid surface: from atomistic
> simulation to thermodynamics*, *Phys. Chem. Chem. Phys.*, 2023, **25**,
> 31352–31362. <https://doi.org/10.1039/D3CP03955B>

The method was originally introduced as *Fields-supported MultiPhase kinetic
Monte Carlo (FsMP/kMC)*.

## Requirements

- A C++ compiler (clang++ is recommended).
- A numerical forcefield (potential) file. See [Forcefields](#forcefields).
- Python 3 (optional) for the post-processing scripts in `xyz_modification/`.

## Building and running

The program is built per configuration: each file in `configs/` includes the
core `program_body.cpp` and defines the parameters of a single simulation.
Use one of the provided configurations (for example `terephthalic_acid_chain.cpp`)
as a template, or write your own.

```bash
cd configs
clang++ -O3 terephthalic_acid_chain.cpp -o program.o
./program.o
```

## Forcefields

The intermolecular interaction is supplied as a precalculated *numerical
potential*. Ready-to-use potentials in the compact binary format (v2) are read by
the program directly; download and unpack them into the `forcefields/` folder:

[Download numerical forcefields (binary, v2)](https://1drv.ms/f/c/18917b5147a88b6c/IgC8SnvBaZORTYWhHGeVLxWQAUzqPePWJhDM3ah1dJotJos?e=5CZNiR)

The original ASCII grids of the DFT potentials are kept in a
[separate folder](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS).
They are only needed to repack a potential yourself, for example with different
folding or in double precision.

The run time reads only the binary format. To convert an ASCII potential (a
legacy one, or your own) use the bundled tool, then point a configuration's
potential path at the resulting `.bin` file:

```bash
clang++ -O3 tools/pack_forcefield.cpp -o pack
./pack forcefields/NAME.dat forcefields/NAME.bin
```

If the molecule has an n-fold rotational symmetry, pass the period in degrees as a
third argument (120 for a C3 molecule, 180 for C2) to store a single period and
shrink the grid; the tool checks the symmetry against the data before folding.

Add `--float` to store the energies in 32-bit precision. The file is half the
size, and the rounding error in the physically relevant region (about 0.01 J/mol)
is negligible compared to the thermal energy.

## Tests

```bash
./tests/run_tests.sh
```

The suite first checks the ASCII-to-binary converter on a synthetic grid
against an independent reimplementation of the packing rules, then runs the
engine on a small grid committed to the repository and compares the
deterministic initial energy of the TMA HCP crystal with a pinned value. When
the full TMA simple potential is present in `forcefields/`, the same check also
runs against the published reference energy.

## Repository layout

| Path | Description |
| --- | --- |
| `program_body.cpp` | Core simulation loop, included by every config file. |
| `configs/` | Ready-to-compile run configurations (compile a file here to build a run). |
| `includes.h` | Master list of headers pulled into `program_body.cpp`. |
| `energies_and_forces_numerical.h` | Intermolecular potential evaluated from the precalculated numerical grid (interpolation, tail correction, hard-core cutoff). |
| `interpolation.h`, `read_forcefield.h` | Grid interpolation and loading of the binary numerical potential. |
| `fields.h` | Damping field, external field, and the pressure change across the gas-solid interface. |
| `Rosenbluth_iteration.h`, `Metropolis_iteration.h` | Kinetic Monte Carlo (Rosenbluth) and Metropolis moves. |
| `StructureGenerator.h` | Generation of the initial molecular structure and unit cell. |
| `pressure_balance.h` | Mechanical equilibrium and pressure balancing. |
| `Widom_test.h` | Widom insertion check of the chemical potential. |
| `Weighted_averages.h`, `block_error.h`, `bootstrap_error.h` | Time averaging and error estimation. |
| `write_xyz_file.h` | Trajectory and configuration output (XYZ). |
| `random/` | SFMT / Mersenne Twister random number generator (by Agner Fog). |
| `forcefields/` | Numerical potential files (downloaded separately). |
| `tools/` | `pack_forcefield.cpp`: converts an ASCII potential into the compact binary grid the run time reads. |
| `tests/` | Regression tests and their small data grid (`./tests/run_tests.sh`). |
| `xyz_modification/` | Python helper scripts for post-processing XYZ trajectories. |

## Status

This is a research code under active cleanup. It reproduces the published
results, but it is not yet packaged for general use. The build is per
configuration (no build system yet), and parts of the optimization routines are
experimental and may be unstable. A regression test suite lives in `tests/`
(see [Tests](#tests)). Improving and documenting the code is ongoing.

## License

Released under the GNU General Public License v3.0. See [LICENSE](LICENSE).

This repository is a fork and continuation of the FSMP-kMC code. The code was
written by S. S. Akimenko and V. A. Gorbunov, under the scientific supervision of
E. A. Ustinov. It is now maintained and further developed by Sergey S. Akimenko.
See the commit history for the changes made in this fork.

## In memory

In memory of Eugene A. Ustinov (1948–2024), whose scientific guidance shaped this
work.
