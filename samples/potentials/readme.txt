The demonstration potential
===========================

TMA_simple_2020_coarse_demo.v2.bin (1 MB) is the only numerical potential
small enough to ship with the program. It is what makes the quickstart work
with no downloads:

    .\fsmp.exe configs\tma_quickstart_demo.txt   (or ./fsmp.out on Linux/macOS)

and it is the potential attached to the ready-to-open Studio project
samples/projects/TMA_quickstart.

What it is
----------
The simplified 2020 trimesic acid model (a single Lennard-Jones center plus
six point charges, see samples/models/TMA_simplified_2020.site), sampled on
every 4th point of the published grid: 281 distances x 31 x 31 angles,
dr = 0.08 A, da = 4 deg, 32-bit values, folded to one 120 deg period.

It is a COARSE grid, kept for demonstration and for the regression suite.
The crystal is stable on it and the numbers land close to the published
ones (initial HCP energy -61.68 kJ/mol against -62.86 on the full grid, a
2 % difference that comes entirely from the coarse sampling), but it is not
a potential to publish numbers from.

For real calculations
---------------------
Download the full potentials into the forcefields/ folder next to this one
and point your configuration or project at them; they are published as a
citable dataset, https://doi.org/10.5281/zenodo.21959125, and the same link
is in forcefields/readme.txt. Or build your own potential: the Studio generates
one from a molecule model on the "Create potential" tab, and
tools/pack_forcefield converts an existing ASCII grid.

Regenerating this file
----------------------
    python3 tests/make_test_grid.py forcefields/TMA_simple_2020.v2.bin \
            samples/potentials/TMA_simple_2020_coarse_demo.v2.bin 4
