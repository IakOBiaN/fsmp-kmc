# Roadmap

Where FSMP-kMC is going, and where it is deliberately not going. Plans change;
this file says what is intended today, not what is promised.

## In flight

- **A permanent home for the potentials.** The precalculated grids are large
  (65 MB to 2 GB each) and currently live in personal cloud folders. They are
  moving to a Zenodo record with a DOI, so they can be cited, verified and
  linked from a paper without depending on one person's account.
- **A citable version of the code.** A Zenodo deposit for the software
  itself, and a submission to the
  [Journal of Open Source Software](https://joss.theoj.org/), which is where
  a piece of research software of this size belongs.

## Considered, not scheduled

- **64-bit iteration counters.** Monte Carlo steps times particles currently
  has to fit in a 32-bit integer, which caps a run at a few million steps for
  a few hundred molecules. The engine refuses such a run instead of
  overflowing, but the ceiling itself should go.
- **Better pair potentials without leaving the app.** The built-in generator
  uses MMFF94, which is honest but underbinds. A neural-network backend
  (AIMNet2 or ANI-2x) would sit behind the same interface as an optional
  install, never bundled: adding a gigabyte of PyTorch to a demonstration app
  would be the wrong trade.
- **More worked examples** reproducing published systems, in the same
  open-and-run form as the bundled quickstart.

## Not planned

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
