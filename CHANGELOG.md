# Changelog

Notable changes per release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [semantic versioning](https://semver.org/spec/v2.0.0.html). What a
major version protects is the data a user accumulates: the binary potential
format and the project file. A key in the parameter file can change in a
minor release, but only when the program refuses the old form with a message
saying what to write instead, and only with an entry here.

Every release is published with ready-made archives for Windows, Linux and
macOS at <https://github.com/IakOBiaN/fsmp-kmc/releases>.

## Unreleased

### Added

- `nnpot`, a standalone command-line tool that computes a pair potential with
  a machine-learned potential (AIMNet2) instead of a classical force field,
  and writes the same v2 binary the engine already reads. The three
  descriptions this project can now put on the trimesic acid dimer place its
  hydrogen bond 38 kJ/mol apart, and that spread propagates into every energy
  the method reports. The grid follows from the molecule and the model alone,
  with nothing calibrated in between. A grid that is interrupted, by Ctrl-C or
  by a reboot, is finished rather than started over on the next run. The tool
  is installed separately and is never part of the Studio or of the release
  bundles, which keep their current size. See `nnpot/README.md`.
- The unit cell optimizer writes the cell it converged on as a `.cell` file,
  next to its xyz animation and named after it. A command-line optimization
  and a production run are two steps of one workflow now: optimize once, then
  point a configuration's `structure` key at the file that came out. The
  Studio has always had that loop, the command line had no way to keep the
  result except by copying numbers out of the log.
- Molecule models for isophthalic and phthalic acid in `samples/models/`,
  next to trimesic and terephthalic acid.
- The published dataset now carries a second set of grids for the same four
  acids, computed with `nnpot` and reaching out to 30 Å, next to the classical
  set it already held. The DOI in the readme resolves to the new version on
  its own, so nothing has to be relinked.

### Changed

- A production run started from the Studio uses the cell on the Unit cell tab
  as it is. It used to hand that cell to the optimizer first, so a deliberately
  imperfect cell, the kind you build to watch a particular behaviour, was
  quietly replaced by the nearest minimum before the Monte Carlo began.
  Optimizing is what the button on tab 4 is for, and it stays there.
- The initial structure now comes from a `.cell` file: a configuration's
  `structure` key holds the path to one, and `samples/cells/` carries the
  sixteen reference cells. The Studio has been reading and writing that same
  format all along, so a cell built there runs from the command line
  unchanged, and a reference cell can be opened, edited and saved back. Until
  now every one of those cells was also written out as C++ inside
  `StructureGenerator.h`, which is where a configuration looked them up by
  name; that second copy is gone, and with it the risk of the two disagreeing.
  **Parameter files that name a structure have to be updated**: replace
  `structure = TMA_HCP_simple_2020` with
  `structure = samples/cells/TMA_HCP_simple_2020.cell`. A name that no longer
  resolves stops the run at once with a message saying exactly that. Nothing
  changes for the Studio, which has always worked in cell files.
  `structure = calculate`, which optimizes the rough cell given by the
  `unit_cell` key, is untouched. One consequence is worth naming: a cell file
  keeps every molecule inside the cell, while the header sometimes placed one
  a whole cell to the left of it. The lattice is identical either way, but the
  two free edges of the starting block are then built from different
  molecules. That is an edge effect and falls off as 1/`uc_in_x`: it reaches
  0.5 kJ/mol per molecule at twenty cells across for the worst of the sixteen
  cells, and 0.0005 for the configurations as they ship. The first Monte Carlo
  steps relax the edge in any case.
- Every reference unit cell in `samples/cells/` was optimized again on the
  corrected potential. Only the honeycomb of the simplified model was already
  at its minimum; the rest had been left behind by an earlier optimizer and
  gain between 0.02 and 2.4 kJ/mol per molecule, the vertical chain of
  isophthalic acid most of all.

### Fixed

- The tail correction had the sign of its slope term flipped, from the commit
  that first wrote it (`c29a8fc`, April 2023) until now. A shifted-force
  truncation subtracts the slope at the cutoff, `u(r) - u(rc) - u'(rc)(r - rc)`,
  so that the force vanishes there. The code added that term instead, which
  doubled the discontinuity of the force rather than removing it and tilted
  every pair energy by a ramp of the wrong sign, growing linearly with the
  distance inside the cutoff. The ramp reaches 3.8 kJ/mol on a single pair of
  the trimesic acid grid, and it accumulates over close neighbours, so it
  favoured dense packings: lattice energies were too deep by 0.2 to 11 kJ/mol
  per molecule, the more tightly packed the structure the more so. Every
  energy the engine reports moves.
- The tail correction read the potential at the nearest angular grid node
  while the energy itself was interpolated between nodes, and the slope term
  multiplied that mismatch by (r - cutoff)/dr, which is of order a thousand
  for close pairs. A pair energy therefore jumped by up to 1.9 kJ/mol whenever
  an orientation crossed the midpoint between two grid nodes: the potential
  was discontinuous in orientation, and the unit-cell optimizer could settle
  on one of those steps instead of on a geometry. Both are now read the same
  way. Energies move by a few thousandths of a kJ/mol on the demonstration
  grid and by up to 0.4 on the published ones; the pinned test energies were
  updated to match.
- The readmes described the published potentials as DFT potentials. They are
  not: the interactions come from the DREIDING force field with its explicit
  hydrogen-bond term, and the quantum calculation behind them produced the
  partial charges and nothing else. The file names always recorded this, `q`
  for the charge calculation and `Dhb` for the hydrogen-bond distance, and the
  prose now agrees with them.
- The documented Windows command line now works in PowerShell too:
  `.\fsmp.exe configs\...` instead of a bare `fsmp.exe configs\...`.
  PowerShell never runs a program from the current folder without the
  prefix, so the quickstart in the release readmes failed for anyone who
  does not use cmd.exe. The same form also survives a machine where
  `NoDefaultCurrentDirectoryInExePath` is set, which makes cmd.exe behave
  the same way.

## v1.1.0 (2026-08-20)

### Added

- A quickstart that runs out of the box: a coarse demonstration potential
  ships in `samples/potentials/`, together with the `TMA_quickstart` project
  and the `configs/tma_quickstart_demo.txt` parameter file. Nothing has to be
  downloaded for a first run, and the regression suite runs the shipped
  configuration as a test.
- **Open the demonstration** on the Studio start page and in the File menu:
  one click installs that project into the user's documents folder and opens
  it on the Run tab. The installed copy carries its own potential, so it keeps
  working when it is moved.
- Every release archive now carries `LICENSE`, `THIRD_PARTY_NOTICES.md` and a
  `licenses/` folder with the terms of everything bundled with it (Qt under
  the LGPL, the Python runtime, NumPy, RDKit, Pillow, PyYAML, the PyInstaller
  bootloader, and SFMT inside the engine). The bundle is not packaged unless
  they are all present.
- A `-cli` archive per platform: the engine, the converter and the example
  data without the Studio, a few megabytes instead of a hundred.
- `BUILD_INFO.txt` in every archive, naming the version, the platform, the
  commit and the frozen library versions it was built from.
- `CITATION.cff`, so the repository offers a citation, and this changelog.
- The full potentials are published as a citable Zenodo dataset
  ([10.5281/zenodo.21959125](https://doi.org/10.5281/zenodo.21959125),
  CC BY 4.0) instead of a personal cloud folder, with a README documenting
  every grid and the binary format.

### Changed

- The README opens with what the program produces, the animation of a real
  run, a quickstart with the numbers to expect, and how the method differs
  from lattice kMC codes.
- The Studio's Run tab keeps the relaxation length from exceeding the run
  length.
- `build.cmd` finds a Python for the test and bundle targets: the Studio
  environment, then `python` on PATH, then the `py` launcher, and it says what
  to install when there is none instead of failing with `python3: not found`.
- One line-ending convention, enforced by `.gitattributes` rather than by each
  machine's git configuration: LF in the repository, and CRLF for batch files,
  which `cmd.exe` cannot run otherwise.

### Fixed

- The engine refuses parameters it cannot run with, naming the file, the line
  and the key: `nan` and `inf`, integers that do not fit in 32 bits, a
  temperature that is not positive, a zero loop step over a non-empty range,
  `free_space` outside `[0, 0.5)`, a relaxation longer than the run, and a
  Monte Carlo step count that would overflow the iteration counter. Before
  this, such a file could run to completion and write `nan` results with a
  zero exit code.

## v1.0.2 (2026-07-21)

### Fixed

- The Windows Studio executable no longer trips antivirus heuristics: the
  release build compiles the PyInstaller bootloader from source instead of
  taking the prebuilt one, which Windows Defender had quarantined as a false
  positive.

## v1.0.1 (2026-07-21)

### Added

- A screenshot tour in the README, version resources and icons on the Windows
  binaries, and `SHA256SUMS.txt` with every release.

### Fixed

- The trajectory viewer survives switching between runs: a stale index no
  longer reaches the newly selected run.
- The macOS instructions match Sequoia (System Settings, Privacy & Security,
  Open Anyway).

## v1.0.0 (2026-07-20)

The whole workflow, documented and shipped.

### Added

- Atomistic models generate their own potentials: MMFF94 through RDKit, with
  geometry optimization in the molecule editor.
- Editing tools for the molecule and the unit cell: center, rotate, aim an
  atom along +x, select and rotate a molecule by a handle.
- Ready-to-open example projects and reference data in `samples/`.
- The Run tab saves its settings on demand and flags a potential file that is
  attached but missing.
- `Studio --selftest`, run on every assembled bundle by the release workflow,
  and a GUI test suite in CI.

### Changed

- MMFF slab generation is about three times faster.
- Cancelling a generation resets the page instead of leaving it stuck.

## v0.5.0 (2026-07-17)

### Added

- Release bundles with FSMP-kMC Studio and the engine for Windows, Linux and
  macOS, assembled by one recipe shared with the local `make bundle`.
- `build.cmd`: the whole Windows toolchain behind one command.
- One version source (`version.h`) for the engine, the Studio, the bundle
  names and the tag guard.

### Changed

- The GUI drives the engine natively on all three platforms; the WSL launcher
  is gone.
- SFMT has a portable branch for non-x86 with the same bit stream.

## v0.4.0 (2026-07-16)

### Added

- The Run tab: detached runs with live progress, statistics plots and a
  trajectory viewer.
- The simulation cell tab, cell files and a folder of reference cells.
- Unit cell optimization driven from the Studio, with a live animation.
- Optimizer stage 0, a global scan of the cell scale, which escapes the capped
  repulsion plateau of a badly scaled starting cell.

### Changed

- The damping field follows the temperature loop; the xyz output separates the
  external field into its own column.

## v0.3.0 (2026-07-13)

### Added

- A stabilization mask that keeps porous polymorphs (chicken-wire, flower)
  from decaying.
- A configurable reference area for the chemical potential (`sigma_mode`) and
  output names that never overwrite an existing file.

### Fixed

- Rotational folding of a numerically asymmetric potential averages the
  symmetric periods instead of keeping one arbitrary period.

## v0.2.0 (2026-07-13)

The release that made the program usable without a compiler.

### Added

- All simulation parameters are read from a text file at run time, so one
  binary serves every system; `configs/` doubles as the annotated template.
- The compact binary potential format (v2) with `tools/pack_forcefield.cpp`,
  optional 32-bit storage and folding by rotational symmetry: 65 GB of ASCII
  grids became 5.2 GB.
- A regression suite, a Makefile, GitHub Actions CI, and a release workflow
  with prebuilt Windows binaries.
- The project logo.

### Changed

- The unit cell optimizer was rewritten: adaptive random descent with
  step adaptation, dead degrees of freedom removed and a convergence report.

### Fixed

- Fixed-size buffers are guarded against overflow, all compiler warnings are
  gone, and a dormant `E_INF` restore bug in the structure generator is fixed.

## v0.1.0 (tagged 2026-07-17)

The inherited baseline: the code as it stood when this fork began, before the
cleanup. It reproduces the published results and is kept as a reference point.
