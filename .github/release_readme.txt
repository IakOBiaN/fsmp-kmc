FSMP-kMC release bundle
=======================

FSMP-kMC Studio   the desktop app - start it and work:
                    Windows:  FSMP-kMC Studio.exe
                    Linux:    ./FSMP-kMC\ Studio
                    macOS:    FSMP-kMC Studio.app (unsigned build: the first
                              time use right-click -> Open)
fsmp.exe / fsmp   the simulation engine, also usable from the command line
pack.exe / pack   converts an ASCII potential to the binary format (optional)
_internal/        the Studio runtime (do not edit)
configs/          example parameter files for command-line runs
models/           molecule models: atomistic .xyz and site .site files
cells/            reference unit cells (.cell) to open in the Studio
forcefields/      put the downloaded potentials here (see readme.txt inside)

Quick start
1. Download the potentials you need (links in forcefields/readme.txt) and
   unpack them into the forcefields folder.
2. Start the Studio and create a project; the bundled models and cells can
   be opened straight from the models and cells folders.

Command line, without the Studio:

    fsmp.exe configs\tma_acid_hcp.txt      (Windows)
    ./fsmp configs/tma_acid_hcp.txt        (Linux, macOS)

Output files are written to the current folder.

Project page, sources and documentation:
  https://github.com/IakOBiaN/fsmp-kmc
