# Roadmap

Where FSMP-kMC is going, and where it is deliberately not going. Plans change;
this file says what is intended today, not what is promised.

## In flight

- **A citable version of the code.** The potentials already have their own
  Zenodo record, and the DOI carried by the badge, the readme and
  `CITATION.cff` resolves to its latest version on its own. The program has
  no such handle. It needs a Zenodo deposit tied to a release tag, and then a
  submission to the
  [Journal of Open Source Software](https://joss.theoj.org/), which is where
  a piece of research software of this size belongs. Everything JOSS asks of
  the repository itself is in place: the license, the CI, the two test
  suites, the contribution guide and a quickstart that needs no download.
  The short paper it reviews is what is missing.

## Considered, not scheduled

- **The surface itself.** The substrate is a plane and nothing more: molecules
  move in it, and no term in the energy knows what they are lying on. For a
  method that reports absolute free energy and entropy that is the first thing
  anyone asks about, and the answer today is that it is not there. The plan is
  a one-body potential V(x, y, alpha), periodic over the substrate cell, added
  through the same external-field slot the damping field and the stabilization
  mask already use, with a Boltzmann average over that cell correcting the
  ideal-gas state the chemical potential is referred to. First as a few
  Fourier terms with period and amplitude as parameters, to find out whether
  it changes anything, then as a tabulated file if it does. Two things have to
  be settled on the way: the box wraps, so the substrate cell has to divide
  it, which continuous pressure balancing does not respect, and the substrate
  takes a share of the stress that the pressure bookkeeping would have to
  gain. What it buys is registry, a physical reason to prefer one site over
  another where the stabilization mask is only a numerical one.
- **64-bit iteration counters.** Monte Carlo steps times particles currently
  has to fit in a 32-bit integer, which caps a run at a few million steps for
  a few hundred molecules. The engine refuses such a run instead of
  overflowing, but the ceiling itself should go.
- **More worked examples** reproducing published systems, in the same
  open-and-run form as the bundled quickstart.

## Not planned

- **The neural potential inside the app.** `nnpot` computes pair potentials
  with AIMNet2 from the command line, and it stays there. PyTorch with CUDA
  is four gigabytes against the sixty megabytes of a release bundle, and a
  workbench that installs like that is a worse program for everyone who never
  asked for a neural potential. Nothing is lost by keeping them apart: the
  grid `nnpot` writes is the same v2 binary as every other potential, and the
  Studio attaches it the same way. The computation is what stays outside.
- **Oblique unit cells.** The elongated two-phase cell and its periodicity
  along y require a rectangular box; tiling an oblique cell leaves a seam
  unless the tilt is commensurate. An oblique lattice is properly expressed as
  a rectangular supercell with a larger molecular basis, which the existing
  chained-polar parameterization already handles.
- **More than one molecular species per project.** It would touch every layer
  from the potential format upwards, and no planned study needs it.
- **Chemistry.** Molecules are rigid bodies by construction: no bonds break,
  form or bend. That assumption is what makes a precalculated pair potential
  possible, and it is not going away.
- **Becoming a general-purpose MD or MC package.** There are good ones
  already.
