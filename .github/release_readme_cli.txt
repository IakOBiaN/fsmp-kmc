FSMP-kMC command-line release
=============================

The simulation engine on its own, without the FSMP-kMC Studio desktop app:
a few megabytes instead of a hundred, nothing to install, no Python, no Qt.
For the Studio, download the full bundle for your platform from the same
release.

fsmp.exe / fsmp   the simulation engine; give it a parameter file
pack.exe / pack   converts an ASCII potential to the binary format (optional)
configs/          example parameter files, each documenting every key
samples/          example data:
                    models/     molecule models (atomistic .xyz and site .site)
                    cells/      reference unit cells (.cell)
                    potentials/ the small demonstration potential the
                                quickstart runs on (see the readme inside)
                    projects/   Studio projects; useful here only as a record
                                of the settings behind each example
forcefields/      put the downloaded production potentials here (see the
                  readme.txt inside)
LICENSE           FSMP-kMC is free software under the GNU GPL v3.0
THIRD_PARTY_NOTICES.md
BUILD_INFO.txt    what this archive was built from

Quick start, with nothing to download

    fsmp.exe configs\tma_quickstart_demo.txt      (Windows)
    ./fsmp configs/tma_quickstart_demo.txt        (Linux, macOS)

Run it from the folder you unpacked, not from inside configs/: paths in a
parameter file are relative to the working directory, and the run writes its
output files there too. Half a minute later you have an xyz trajectory, the
optimized unit cell and a statistics table for a trimesic acid monolayer at
300 K, computed on the bundled demonstration grid.

That grid is deliberately coarse, for a fast first run. For real numbers,
download a full potential (links in forcefields/readme.txt), unpack it into
the forcefields folder and point a configuration's `potential` key at it.

The binaries are static: no runtime libraries are needed. Windows
SmartScreen may still warn about an unrecognized publisher on the first
launch of a downloaded file (More info -> Run anyway); on macOS 15+ allow it
under System Settings -> Privacy & Security -> Open Anyway, or unpack with
`tar -xzf` in Terminal, which never sets the quarantine flag.

Project page, sources and documentation:
  https://github.com/IakOBiaN/fsmp-kmc
