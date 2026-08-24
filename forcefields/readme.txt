Unpack the potential archives into this folder.

You do not need any of them to try the program: the quickstart
(configs/tma_quickstart_demo.txt and the Studio project
samples/projects/TMA_quickstart) runs on the small demonstration grid in
samples/potentials/. Download these for real calculations.

Ready-to-use binary potentials (format v2, read by the program directly),
published as a citable dataset under CC BY 4.0:
https://doi.org/10.5281/zenodo.21959125

  TMA_simple_2020.v2.bin                       66 MB   trimesic acid, simplified model
  TMA_q(B3LYP+PBE)_...v2.bin                  302 MB   trimesic acid, atomistic
  TPA_qPBE_crystal_...v2.bin                  548 MB   terephthalic acid
  IPA_qPBE_crystal_...v2.bin                  2.2 GB   isophthalic acid
  PA_qPBE_crystal_...v2.bin                   2.2 GB   phthalic acid

The atomistic grids were computed with the DREIDING force field and its
explicit hydrogen-bond term; the quantum calculation behind them produced the
partial charges and nothing else. A name records both: q<method> is the
charge calculation, Dhb the hydrogen-bond distance.

The same record also holds a second set for the same four acids, computed
with the AIMNet2 machine-learned potential and reaching out to 30 A:

  TMA_aimnet2_...v2.bin                        66 MB   trimesic acid
  TPA_aimnet2_...v2.bin                       161 MB   terephthalic acid
  IPA_aimnet2_...v2.bin                       639 MB   isophthalic acid
  PA_aimnet2_...v2.bin                        639 MB   phthalic acid

They are a different description of the same molecules, not a correction of
the first set, and the two disagree by 15 kJ/mol on the trimesic acid dimer.
They were produced with nnpot/, see nnpot/README.md.

Keep the file names as they are: the example configurations and the sample
Studio projects point at them.

Original ASCII grids of these potentials (only needed to repack a potential
yourself, for example with different folding or in double precision; see
tools/pack_forcefield):
https://1drv.ms/f/s!AmyLqEdRe5EYgdkXdo7VUsFQxyMmng?e=6Vi3NS
