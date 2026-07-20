FSMP-kMC release bundle
=======================

FSMP-kMC Studio   the desktop app - start it and work:
                    Windows:  FSMP-kMC Studio.exe (SmartScreen may warn on
                              the first launch: More info -> Run anyway.
                              Windows Defender may also quarantine the exe
                              as a false positive - it is the PyInstaller
                              packer, not malware; restore it from Windows
                              Security -> Protection history and exclude the
                              folder before re-extracting)
                    Linux:    ./FSMP-kMC\ Studio (needs glibc 2.35+,
                              i.e. Ubuntu 22.04 or newer)
                    macOS:    FSMP-kMC Studio.app (unsigned build; macOS 15+
                              blocks the first launch: allow it in System
                              Settings -> Privacy & Security -> Open Anyway;
                              on older macOS right-click -> Open is enough.
                              Unpacking with `tar -xzf` in Terminal avoids
                              the prompt entirely)
fsmp.exe / fsmp   the simulation engine, also usable from the command line
pack.exe / pack   converts an ASCII potential to the binary format (optional)
_internal/        the Studio runtime (do not edit)
configs/          example parameter files for command-line runs
samples/          example data to explore:
                    models/    molecule models (atomistic .xyz and site .site)
                    cells/     reference unit cells (.cell) to open in the Studio
                    projects/  ready-to-open Studio projects (open one and run)
forcefields/      put the downloaded potentials here (see readme.txt inside)

Quick start
1. Download the potentials you need (links in forcefields/readme.txt) and
   unpack them into the forcefields folder.
2. Start the Studio and either open a ready-made project from samples/projects,
   or create your own from the bundled samples/models and samples/cells.

Command line, without the Studio:

    fsmp.exe configs\tma_acid_hcp.txt      (Windows)
    ./fsmp configs/tma_acid_hcp.txt        (Linux, macOS)

Output files are written to the current folder.

Project page, sources and documentation:
  https://github.com/IakOBiaN/fsmp-kmc
