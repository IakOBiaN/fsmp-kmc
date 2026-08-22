# Neural-network pair potentials

`nnpot` builds an FSMP-kMC pair potential with a machine-learned potential
instead of a classical force field. It sweeps two rigid copies of a molecule
over the grid of centre-to-centre distances and orientations, asks a neural
network for the energy of every arrangement, and writes the same binary `v2`
file the engine already reads. The result is attached to a project like any
other potential.

The tool lives outside FSMP-kMC Studio on purpose. PyTorch with CUDA is about
four gigabytes, and the Studio bundle is sixty megabytes and stays that way.
Nothing here is imported by the Studio, by the engine, or by the release
bundle: this is a separate program that produces a file.

## Why bother

The trimesic acid dimer at its hydrogen-bonded minimum, scored by the two
generators this repository ships, from the same molecule and the same
geometry:

| Generator | Well depth | Distance |
| --- | --- | --- |
| AIMNet2, this tool | -63.2 kJ/mol | 9.72 Å |
| MMFF94, the Studio generator | -40.0 kJ/mol | 9.96 Å |

Both find the same arrangement, the double hydrogen bond between two
carboxyls, and they disagree about its strength by 23 kJ/mol. That gap is not
cosmetic: the well depth feeds straight into the cohesion energy, the chemical
potential and the transition temperatures the method exists to produce.

Nothing here is calibrated against anything. The grid is what the model says
about the molecule you hand it, and that is the property worth having: the
same molecule and the same model give the same numbers on any machine, with
no step in between that someone has to trust.

## Installation

The environment is separate from `gui/.venv` deliberately, so that PyTorch
can never be picked up by the bundle builder and shipped to users.

```bash
py -3.11 -m venv nnpot/.venv
nnpot/.venv/Scripts/python -m pip install -r nnpot/requirements.txt
```

On Linux and macOS the interpreter is `nnpot/.venv/bin/python` instead.

The versions in `requirements.txt` are pinned for a reason. AIMNet2 needs
`torch-cluster`, a compiled extension, and prebuilt wheels for it exist only
for some combinations of PyTorch version, CUDA version, Python version and
platform. PyTorch 2.8 with CUDA 12.8 is the newest combination that has one
for Windows and Python 3.11; with a newer PyTorch the extension would have to
be built from source. On Linux, or with a compiler at hand, other
combinations work equally well.

The model weights are downloaded on first use and cached inside the
environment.

The tool imports the `v2` file writer from `gui/fsmp_gui`, so it runs from a
checkout of this repository rather than as a lone script. That is deliberate:
the binary format and the folding rules have exactly one implementation,
shared with the Studio and matching `tools/pack_forcefield.cpp`.

## Use

Relax the molecule with the model that will compute the potential, so that the
geometry and the interaction come from the same place and the result depends
on nothing but the molecule you started from. It is kept in its plane, because
the engine's molecules are planar rigid bodies:

```bash
python -m nnpot optimize --molecule samples/models/trimesic_acid.xyz --out tma_aimnet2.xyz --symmetrize
```

`--symmetrize` makes the rotational symmetry exact. It is worth using: an
optimizer leaves a threefold molecule a few thousandths of an ångström short
of exact C3, and folding the grid to one 120 degree period is only legitimate
when the periods really are identical.

Then compute the grid:

```bash
python -m nnpot grid --molecule tma_aimnet2.xyz --out TMA_aimnet2.v2.bin
```

The defaults produce a production grid: distances from the molecule's own
radius out to 30 Å in steps of 0.02 Å, angles every degree, folded to the
detected rotational symmetry, stored as float32. The options that matter:

| Option | Meaning |
| --- | --- |
| `--r-min`, `--r-max`, `--dr` | distance range and step, ångström |
| `--da` | angular step, degrees |
| `--fold` | stored angular period; defaults to 360 over the symmetry |
| `--wall` | closest atom-atom contact the model is trusted with, ångström |
| `--batch` | arrangements per call; defaults to a safe size for the molecule |
| `--double` | store float64 instead of float32 |
| `--symmetrize` | make the molecule exactly symmetric first |
| `--charge` | net molecular charge |
| `--model` | model name, default `aimnet2` |
| `--coulomb` | `simple`, `dsf` or `ewald` |

The geometry that went into the grid is part of the result. A parameter file
that uses the potential should point `molecule_model` at the same relaxed and
symmetrized `.xyz`, not at the one it was derived from, because the molecular
area, the hard core and every drawing follow from it.

Interrupting with Ctrl-C finishes the current distance, removes the partial
file and exits. There is no resume: a full grid takes about half an hour, and
adding one would mean a second implementation of the file format.

## Comparing with another grid

```bash
python -m nnpot compare --grid TMA_aimnet2.v2.bin --reference another.v2.bin
```

This compares two grids point by point on whatever they have in common: the
overlapping distances, the coarser of the two steps, and the shared folded
period. Two grids computed from different geometries usually disagree about
where zero degrees is, and one may be the mirror image of the other, so the
comparison searches for the turn and the reflection that bring them into best
agreement and reports which it used. It then gives the root mean square
difference, the systematic bias and the worst single point, counting only
arrangements below `--threshold`, because otherwise the repulsive wall would
dominate every number.

What this is for is worth stating plainly. Comparing against a published grid
tells you how two descriptions differ. It does not tell you which one is
right. Any hand-computed grid carries its own choices, a functional, a basis
set, a charge model, an empirical correction, and none of them arrive with a
guarantee. Read the output as information about the difference, never as an
error to be minimised, and do not reach for a different geometry or different
settings here until the numbers line up with someone else's. That would be
fitting, and it would quietly destroy the one property that makes this tool
worth having: that the grid follows from the molecule and the model alone.

## What the tool does beyond calling the model

**Symmetry.** A grid of one-degree angles over a full turn is 130 321
arrangements per distance. Two exact symmetries cut that to 7 260, a factor
of 18: the molecule's own rotational symmetry, and the fact that swapping
which molecule sits at the origin cannot change the energy of a pair. Both
are verified numerically by the test suite rather than assumed.

**A wall at short range.** Neural potentials are trained on sensible
geometries. Push two atoms into each other and the model can return a deep
well where physics demands a wall, and a spurious well is not a small error
in kinetic Monte Carlo: the system collapses into it. Arrangements whose
closest atom-atom contact is below `--wall` are given the repulsive cap
without ever being shown to the model.

This is not a hypothetical. Sweeping the trimesic acid dimer over the
arrangements the wall rejects and asking the model anyway, about a quarter of
a percent of them come back attractive, the deepest at -29 kJ/mol with two
atoms 0.445 Å apart. That is a rare event and a modest energy next to the real
minimum, but a well that does not exist is a trap a long Monte Carlo run will
eventually find and settle into.

**Long range.** A pair potential for this method has to stay meaningful out
to 30 Å, where the interaction is only a few joules per mole but the virial
and the pressures depend on it. Most machine-learned potentials are strictly
short-ranged and would leave a hard zero beyond about 14 Å. AIMNet2 is not:
it predicts partial charges and adds an unscreened Coulomb term over all
pairs, so its tail is real, and the tool checks that the model says so rather
than assuming it. For a model that does declare a cutoff, the tool adds the
Coulomb interaction of the pairs beyond that cutoff, switched on smoothly, so
that the two descriptions join without a seam and without a double count.

**Memory.** The sizes of the tensors inside the model depend on the geometry,
so PyTorch's caching allocator grows steadily instead of reusing blocks. On
Windows it then spills into system memory and keeps working, silently, about
twenty-five times slower. The tool watches the reserved memory and releases
the cache when it passes half the card. The coarse grid below takes 28
seconds with that guard and 10.2 minutes without it.

## Cost

Measured on a laptop RTX 3080, trimesic acid, 21 atoms:

| Grid | Arrangements | Time |
| --- | --- | --- |
| 281 by 31 by 31, step 0.08 Å, 4 degree angles | 126 thousand | 28 s |
| 1121 by 121 by 121, step 0.02 Å, 1 degree angles | 8.1 million | 37 min |

The second is the size of the published production grids. Throughput runs
between 3 500 and 4 400 arrangements per second depending on the grid.

## Tests

```bash
python -m unittest discover -s nnpot/tests
```

Twelve of the tests need only NumPy and RDKit, so they run in continuous
integration alongside the Studio suite and PyTorch never enters it. They
check the machinery against the Studio's own MMFF94 generator, and the
decisive one compares the files byte for byte: the same grid computed through
this tool and through the Studio has to be identical. That is what makes the
symmetry shortcuts and the file writing trustworthy.

They also cover the comparison: the same grid computed from a molecule
turned by 30 degrees has to be recognised as the same grid, turned.

Two more tests exercise the model itself and are skipped when AIMNet2 is not
installed. In the environment above the whole suite is fourteen tests.

## Limits

Molecules are planar rigid bodies, which is the engine's model rather than a
limitation of this tool alone. The relaxation keeps the molecule in its
plane, so a species whose real minimum is twisted is being approximated.

AIMNet2 covers H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br and I. The tool
refuses anything outside the elements it can map.

The accuracy of the result is the accuracy of the model. The agreement above
is one molecule against one reference; a new system deserves its own check
against whatever higher-level data exists for it.
