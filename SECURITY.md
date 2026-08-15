# Security policy

## What this program is exposed to

FSMP-kMC is an offline simulation program. It opens no sockets, contacts no
server and needs no elevated rights. Its attack surface is the files it
reads: parameter files, molecule models (`.xyz`, `.site`), unit cells
(`.cell`), Studio project files and, above all, binary potential grids, which
routinely arrive from someone else as a multi-gigabyte download.

Treat a potential file from an untrusted source the way you would treat any
other untrusted binary. A crafted file that makes the engine read out of
bounds is a bug worth reporting; the parsers check what they can (the grid
header is validated against the file size, the parameter file is rejected on
anything it cannot run with), but they were written to catch honest mistakes,
not deliberate attacks.

## Supported versions

Fixes go into the next release from the current line. Only the
[latest release](https://github.com/IakOBiaN/fsmp-kmc/releases/latest) is
supported; older bundles are kept for reproducibility, not maintained.

## Reporting

Report a suspected vulnerability privately through
[GitHub's security advisories](https://github.com/IakOBiaN/fsmp-kmc/security/advisories/new),
not in a public issue. Please include what you fed the program, what it did,
and the `BUILD_INFO.txt` of the archive you used. Expect an acknowledgement
within a week; this is a one-person project, so a fix may take longer, and
you will be credited in the advisory unless you ask otherwise.

Crashes and wrong numbers that do not involve a hostile input are ordinary
bugs: please open a normal issue instead, where the discussion can be public.
