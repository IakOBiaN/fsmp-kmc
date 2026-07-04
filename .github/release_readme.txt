FSMP-kMC - Windows build
========================

fsmp.exe      the simulation program
pack.exe      converts an ASCII potential to the binary format (optional)
configs/      example parameter files - copy one and edit it for your system
forcefields/  put the downloaded potentials here (see readme.txt inside)

Quick start
1. Download the potentials you need (links in forcefields/readme.txt) and
   unpack them into the forcefields folder.
2. Open the command prompt in this folder and run:

     fsmp.exe configs\tma_acid_hcp.txt

   Output files are written to the current folder.

Project page, sources and documentation:
  https://github.com/IakOBiaN/fsmp-kmc
