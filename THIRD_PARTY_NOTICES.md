# Third-party notices

FSMP-kMC itself is released under the GNU General Public License v3.0; the
full text is in [LICENSE](LICENSE), and it ships in every release archive.
This file lists the third-party components that a release archive carries
alongside it, and the terms they come under.

The `licenses/` folder of a release archive holds the license texts of the
bundled components. Where a component ships its own text, that copy is taken
straight from the version that was frozen into the archive, so the text and
the binary always match.

## In the engine (`fsmp`, `pack`, and every build from source)

**SFMT random number generator** — `random/`
Algorithm by Mutsuo Saito and Makoto Matsumoto (Hiroshima University); the
C++ implementation used here is by Agner Fog, <https://www.agner.org/random/>.
Agner Fog's implementation is under the **GNU General Public License**, and
it contains parts of the original C code published under a **3-clause BSD
license**, which is therefore in effect in addition to the GPL. The GPL text
is `LICENSE` (identical to `random/license.txt` in the source tree); the BSD
notice is reproduced here as required for binary distribution:

```
Copyright (c) 2006, 2007 by Mutsuo Saito, Makoto Matsumoto and Hiroshima University.
Copyright (c) 2008 by Agner Fog.
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    > Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    > Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    > Neither the name of the Hiroshima University nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

The Windows binaries are linked statically against the **MinGW-w64 runtime**
(w64devkit), whose runtime libraries carry permissive licenses with no
attribution requirement for binaries; see
<https://www.mingw-w64.org/licensing/>.

## In the Studio bundle

This section applies to the full bundle only. The `-cli` archive of a
release carries the engine and the example data, nothing from the list
below.

The desktop application is frozen with PyInstaller, so a release archive
contains a Python runtime and the libraries below.

| Component | License | Text |
| --- | --- | --- |
| Qt 6, via **PySide6** and **shiboken6** | GNU Lesser General Public License v3.0 | `licenses/LGPL-3.0.txt` |
| **Python** runtime | Python Software Foundation License 2.0 | <https://docs.python.org/3/license.html> |
| **NumPy** | BSD 3-Clause | `licenses/numpy/` |
| **RDKit** | BSD 3-Clause | `licenses/rdkit/` |
| **Pillow** | MIT-CMU | `licenses/pillow/` |
| **PyYAML** | MIT | `licenses/PyYAML/` |
| **PyInstaller** bootloader | GPL 2.0 with the PyInstaller bootloader exception (a frozen application may carry any license) | `licenses/pyinstaller/` |
| **OpenSSL** (through the Python runtime) | Apache License 2.0 | <https://www.openssl.org/source/license.html> |
| **SQLite** (through the Python runtime) | public domain | <https://www.sqlite.org/copyright.html> |
| **libffi**, **expat**, **zlib**, **bzip2**, **XZ Utils** (through the Python runtime) | MIT, MIT, zlib, BSD-4-clause-like, 0BSD respectively | <https://docs.python.org/3/license.html> |
| **Microsoft Visual C++ runtime** (`VCRUNTIME140*.dll`, `MSVCP140*.dll`, Windows archives only) | Microsoft redistributable terms | <https://visualstudio.microsoft.com/license-terms/> |

### Qt and the LGPL

The Studio uses Qt 6 as a shared library through PySide6, under the LGPL
v3.0. Nothing in Qt is modified. Because FSMP-kMC is itself free software
under the GPL v3.0 and its complete source is published at
<https://github.com/IakOBiaN/fsmp-kmc>, a recipient can rebuild the
application against a different, compatible version of Qt: install the
project's Python package (`pip install ./gui`, which fetches PySide6 from
PyPI, or point it at your own Qt build) and run it from source, or rebuild
the frozen bundle with `python tools/make_bundle.py`. The Qt sources
themselves are available from <https://download.qt.io/>.

The PySide6 wheels also carry a reference to the alternative Qt commercial
license. It does not apply here: these builds use Qt under the LGPL v3.0.
