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
                    models/     molecule models (atomistic .xyz and site .site)
                    cells/      reference unit cells (.cell) to open in the Studio
                    potentials/ the small demonstration potential the
                                quickstart runs on (see the readme inside)
                    projects/   ready-to-open Studio projects; TMA_quickstart
                                works with nothing downloaded
forcefields/      put the downloaded production potentials here (see the
                  readme.txt inside)
LICENSE           FSMP-kMC is free software under the GNU GPL v3.0
THIRD_PARTY_NOTICES.md
licenses/         what this archive carries besides FSMP-kMC itself (Qt,
                  Python, NumPy, RDKit and the rest) and under what terms

Quick start, with nothing to download
1. Start the Studio and press "Open the demonstration". It copies the
   ready-made TMA_quickstart project into your documents folder, opens it
   and lands on the Run tab: press Start there. Half a minute later you have
   a trajectory, live plots and a statistics table for a trimesic acid
   monolayer at 300 K.
   The same run from the command line:
       .\fsmp.exe configs\tma_quickstart_demo.txt    (Windows)
       ./fsmp configs/tma_quickstart_demo.txt        (Linux, macOS)
2. For real numbers, download a full potential (links in
   forcefields/readme.txt), unpack it into the forcefields folder and open
   one of the other sample projects, or build your own model in the Studio:
   the quickstart potential is a coarse demonstration grid.

Command line, without the Studio: every file in configs/ is a documented
parameter file, and the engine takes one as its only argument.

    .\fsmp.exe configs\tma_acid_hcp.txt    (Windows)
    ./fsmp configs/tma_acid_hcp.txt        (Linux, macOS)

(the leading .\ matters: PowerShell will not run a program from the
current folder without it, and cmd.exe accepts it too)

That one is the published trimesic acid setup and needs the full potential
downloaded into forcefields/; configs/tma_quickstart_demo.txt does not.
Output files are written to the current folder.

Project page, sources and documentation:
  https://github.com/IakOBiaN/fsmp-kmc
