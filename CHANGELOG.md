# Changelog

Notable changes per release. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [semantic versioning](https://semver.org/spec/v2.0.0.html): the
parameter file, the binary potential format and the project file are the
interfaces that a major version protects.

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
- Molecule models for isophthalic and phthalic acid in `samples/models/`,
  next to trimesic and terephthalic acid.
- The published dataset now carries a second set of grids for the same four
  acids, computed with `nnpot` and reaching out to 30 Å, next to the classical
  set it already held. The DOI in the readme resolves to the new version on
  its own, so nothing has to be relinked.

### Fixed

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
