# Contributing to FSMP-kMC

Bug reports, questions and patches are all welcome. This is a small
scientific project maintained by one person, so please expect replies in days
rather than hours.

If you are reporting something that does not work, or a result that looks
wrong, go straight to [the issue tracker](https://github.com/IakOBiaN/fsmp-kmc/issues):
there is a template for a plain bug and a separate one for a physics result
that does not look right, because the two are investigated very differently.
A `BUILD_INFO.txt` sits in every release archive; pasting it saves a round of
questions.

## Building and testing

The engine is standard C++ with no dependencies:

```bash
make                     # fsmp.out
make test                # the regression suite (needs Python 3)
```

On Windows, `build` and `build test` do the same through the w64devkit
toolchain without touching `PATH`. The suite runs on the demonstration
potential committed to `samples/potentials/`, so it works right after a
clone; the block that needs a full potential skips itself when the file is
absent.

The Studio is a Python package in `gui/`:

```bash
python3 -m venv gui/.venv
gui/.venv/bin/pip install -e gui
QT_QPA_PLATFORM=offscreen gui/.venv/bin/python -m unittest discover -s gui/tests
```

Both suites run in CI on every push, on Linux, macOS and native Windows.

## What a patch has to keep intact

**The pinned energies.** The suite checks numbers that were validated against
published results: the initial energy of the TMA HCP crystal on the
demonstration grid (-61.7449 kJ/mol) and on the full potential (-62.8605), the
optimizer's final energy (-62.2276), and the shipped quickstart (-61.6841). A
change that moves them is either a bug or a genuine improvement in the
physics; either way it needs an explanation in the pull request, not a quietly
updated constant.

**The warning gate.** CI compiles with `-Wall -Wextra -Werror` under both GCC
and Clang:

```bash
g++ -O2 -Wall -Wextra -Werror -fsyntax-only fsmp.cpp
g++ -O2 -Wall -Wextra -Werror -fsyntax-only -DSFMT_NO_SSE2 fsmp.cpp
g++ -O2 -Wall -Wextra -Werror -fsyntax-only tools/pack_forcefield.cpp
```

**The command line.** The Studio is a front end for the engine, never a
requirement: every workflow must stay available to someone running `fsmp` on a
parameter file over ssh.

**The interfaces.** The parameter file, the binary potential format (v2) and
the Studio project file are what other people's work depends on. Extending
them is fine; changing the meaning of something that already exists is a
major-version matter.

## Things worth knowing before you start

- The engine has no dependencies, and it should stay that way. The Studio
  depends on PySide6, NumPy and RDKit; a new dependency there needs a good
  reason, since it is frozen into a bundle for every platform.
- Do not commit potentials or other large binaries. The one exception is the
  1 MB demonstration grid in `samples/potentials/`, which the tests and the
  quickstart both use.
- Shell scripts and text files use LF endings (`.gitattributes` enforces it).
- Commit messages in this repository are short sentences saying what the
  change does for a user, not what was edited.
- Scratch work belongs in `experiments/`, which is ignored by git.

## Licensing

FSMP-kMC is GPL-3.0. By contributing you agree that your contribution is
released under the same license. There is no contributor licence agreement to
sign.
