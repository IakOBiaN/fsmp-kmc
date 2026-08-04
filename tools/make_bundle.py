#!/usr/bin/env python3
"""Assemble a release bundle: the Studio (PyInstaller), the engine binaries
and the bundled data, laid out exactly like a published release archive.

    python3 tools/make_bundle.py [--name v0.5.0] [--build-engine]

The engine binaries (fsmp.exe/pack.exe on Windows, fsmp/pack elsewhere) are
expected in the repository root: `make windows` builds them on Windows,
--build-engine compiles them here with the release flags. PyInstaller is
looked up on PATH, then in gui/.venv (pip install ./gui[build] provides it).
The bundle folder and its archive land in dist/. The release workflow calls
this same script, so the local and the published layouts cannot drift.
"""

import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
from importlib import metadata
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXE = ".exe" if os.name == "nt" else ""

# Components frozen into the Studio that ship their own license text. It is
# taken from the installed distribution, so the text in the archive always
# describes the version that was actually bundled. Qt is not here: the
# PySide6 wheels carry only a reference to the commercial license, while
# these builds use Qt under the LGPL, whose text is kept in licenses/.
WHEEL_LICENSES = ("numpy", "rdkit", "pillow", "PyYAML", "pyinstaller")


def project_version() -> str:
    match = re.search(r'#define\s+FSMP_VERSION\s+"([^"]+)"',
                      (REPO / "version.h").read_text(encoding="utf-8"))
    return match.group(1) if match else "unknown"


def platform_tag() -> str:
    system = {"Windows": "windows", "Linux": "linux",
              "Darwin": "macos"}.get(platform.system(), sys.platform)
    machine = platform.machine().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return f"{system}-{machine}"


def run(cmd) -> None:
    print("+", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=REPO, check=True)


def find_pyinstaller() -> str:
    found = shutil.which("pyinstaller")
    if found:
        return found
    venv = REPO / "gui" / ".venv"
    for candidate in (venv / "Scripts" / "pyinstaller.exe",
                      venv / "bin" / "pyinstaller"):
        if candidate.is_file():
            return str(candidate)
    sys.exit("pyinstaller not found: pip install ./gui[build]")


def build_engine() -> None:
    static = [] if platform.system() == "Darwin" else ["-static"]
    cxx = os.environ.get("CXX", "g++")
    run([cxx, "-O3", *static, "-Wall", "-Wextra",
         "fsmp.cpp", "-o", "fsmp" + EXE])
    run([cxx, "-O3", *static, "-Wall", "-Wextra",
         Path("tools") / "pack_forcefield.cpp", "-o", "pack" + EXE])


def wheel_license_files(name: str) -> list[tuple[Path, Path]]:
    """(path inside the .dist-info, absolute path) of every license text an
    installed distribution ships."""
    try:
        dist = metadata.distribution(name)
    except metadata.PackageNotFoundError:
        return []
    found = []
    for entry in dist.files or []:
        parts = Path(str(entry)).parts
        info = next((i for i, p in enumerate(parts)
                     if p.endswith(".dist-info")), None)
        if info is None:
            continue                     # package data, not a license text
        rest = parts[info + 1:]
        if not rest:
            continue
        # the folder test comes first: "licenses" itself starts with "license"
        if rest[0] == "licenses" and len(rest) > 1:
            rel = Path(*rest[1:])      # wheels nest their texts one level down
        elif rest[0].upper().startswith(("LICENSE", "COPYING", "NOTICE")):
            rel = Path(*rest)
        else:
            continue
        path = Path(dist.locate_file(entry))
        if path.is_file():
            found.append((rel, path))
    return found


def copy_licenses(bundle: Path) -> None:
    """The project's own license, the third-party notices, and the license
    texts of everything frozen into the archive. Distributing binaries
    without them would breach both the project's GPL and the terms of the
    bundled libraries."""
    shutil.copy2(REPO / "LICENSE", bundle / "LICENSE")
    shutil.copy2(REPO / "THIRD_PARTY_NOTICES.md", bundle)
    dest = bundle / "licenses"
    dest.mkdir(exist_ok=True)
    for text in sorted((REPO / "licenses").glob("*.txt")):
        shutil.copy2(text, dest / text.name)
    for name in WHEEL_LICENSES:
        files = wheel_license_files(name)
        if not files:
            print(f"  ! no license text found for {name}", flush=True)
            continue
        for rel, path in files:
            target = dest / name / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, target)


def studio_name() -> str:
    if platform.system() == "Darwin":
        return "FSMP-kMC Studio.app"
    return "FSMP-kMC Studio" + EXE


def verify(bundle: Path) -> None:
    """Refuse to publish an archive that lost something a user needs: the
    licenses, the app, the engine, or the data the demonstration runs on."""
    required = [
        studio_name(),
        "fsmp" + EXE, "pack" + EXE,
        "README.txt", "LICENSE", "THIRD_PARTY_NOTICES.md",
        "licenses/LGPL-3.0.txt",
        "forcefields/readme.txt",
        "configs/tma_quickstart_demo.txt",
        "samples/models/trimesic_acid.xyz",
        "samples/projects/TMA_quickstart/project.json",
    ]
    missing = [name for name in required if not (bundle / name).exists()]
    # the demonstration is useless without its grid; the name is free to
    # change, the folder may not be empty
    if not list((bundle / "samples" / "potentials").glob("*.v2.bin")):
        missing.append("samples/potentials/*.v2.bin")
    if missing:
        sys.exit("bundle is incomplete, refusing to package:\n  "
                 + "\n  ".join(missing))
    print(f"verified {len(required) + 1} required entries")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default=None,
                        help="version part of the bundle name (e.g. v0.5.0);"
                             " defaults to version.h")
    parser.add_argument("--build-engine", action="store_true",
                        help="compile the engine here with the release flags")
    args = parser.parse_args()
    name = args.name or ("v" + project_version())

    if args.build_engine:
        build_engine()
    engine = REPO / ("fsmp" + EXE)
    pack = REPO / ("pack" + EXE)
    if not (engine.is_file() and pack.is_file()):
        sys.exit(f"engine binaries not found in {REPO}: run `make windows` "
                 "(Windows) or rerun with --build-engine")

    run([find_pyinstaller(), "--noconfirm", "--distpath", "gui/dist",
         "--workpath", "gui/build", "gui/studio.spec"])

    bundle = REPO / "dist" / f"fsmp-kmc-{name}-{platform_tag()}"
    if bundle.exists():
        shutil.rmtree(bundle)
    bundle.mkdir(parents=True)

    if platform.system() == "Darwin":
        # the data folders sit next to the app bundle (engine.app_root)
        shutil.copytree(REPO / "gui" / "dist" / "FSMP-kMC Studio.app",
                        bundle / "FSMP-kMC Studio.app", symlinks=True)
    else:
        shutil.copytree(REPO / "gui" / "dist" / "studio", bundle,
                        dirs_exist_ok=True)
    shutil.copy2(engine, bundle)
    shutil.copy2(pack, bundle)
    # configs are the engine's example parameter files; samples carries the
    # example molecule models, unit cells and ready-to-open Studio projects
    for folder in ("configs", "samples"):
        shutil.copytree(REPO / folder, bundle / folder)
    (bundle / "forcefields").mkdir()
    shutil.copy2(REPO / "forcefields" / "readme.txt", bundle / "forcefields")
    shutil.copy2(REPO / ".github" / "release_readme.txt",
                 bundle / "README.txt")
    copy_licenses(bundle)
    verify(bundle)

    kind = "zip" if os.name == "nt" else "gztar"
    archive = shutil.make_archive(str(bundle), kind,
                                  root_dir=bundle.parent,
                                  base_dir=bundle.name)
    print("bundle: ", bundle)
    print("archive:", archive)


if __name__ == "__main__":
    main()
