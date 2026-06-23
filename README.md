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
potential*. Download the available potentials and unpack them into the
`forcefields/` folder:

[Download numerical forcefields](https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS)

## Repository layout

| Path | Description |
| --- | --- |
| `program_body.cpp` | Core simulation loop, included by every config file. |
| `configs/` | Ready-to-compile run configurations (compile a file here to build a run). |
| `includes.h` | Master list of headers pulled into `program_body.cpp`. |
| `energies_and_forces_*.h` | Alternative implementations of the intermolecular potential (exact, approximate, numerical Dreiding TMA, simple model). |
| `interpolation.h`, `read_forcefield.h` | Grid interpolation and loading of the precalculated numerical potential. |
| `fields_*.h` | Damping and external field definitions (per publication). |
| `Rosenbluth_iteration.h`, `Metropolis_iteration.h` | Kinetic Monte Carlo (Rosenbluth) and Metropolis moves. |
| `StructureGenerator*.h` | Generation of the initial molecular structure and unit cell. |
| `pressure_balance.h` | Mechanical equilibrium and pressure balancing. |
| `Widom_test.h` | Widom insertion check of the chemical potential. |
| `Weighted_averages.h`, `block_error.h`, `bootstrap_error.h` | Time averaging and error estimation. |
| `write_xyz_file.h` | Trajectory and configuration output (XYZ). |
| `random/` | SFMT / Mersenne Twister random number generator (by Agner Fog). |
| `forcefields/` | Numerical potential files (downloaded separately). |
| `xyz_modification/` | Python helper scripts for post-processing XYZ trajectories. |

## Status

This is a research code under active cleanup. It reproduces the published
results, but it is not yet packaged for general use. The build is per
configuration (no build system yet), there is no automated test suite, and parts
of the optimization routines are experimental and may be unstable. Improving and
documenting the code is ongoing.

## License

Released under the GNU General Public License v3.0. See [LICENSE](LICENSE).

This repository is a fork and continuation of the FSMP-kMC code. The code was
written by S. S. Akimenko and V. A. Gorbunov, under the scientific supervision of
E. A. Ustinov. It is now maintained and further developed by Sergey S. Akimenko.
See the commit history for the changes made in this fork.

## In memory

In memory of Eugene A. Ustinov (1948–2024), whose scientific guidance shaped this
work.
