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
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXE = ".exe" if os.name == "nt" else ""


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default="local",
                        help="version part of the bundle name (e.g. v0.5.0)")
    parser.add_argument("--build-engine", action="store_true",
                        help="compile the engine here with the release flags")
    args = parser.parse_args()

    if args.build_engine:
        build_engine()
    engine = REPO / ("fsmp" + EXE)
    pack = REPO / ("pack" + EXE)
    if not (engine.is_file() and pack.is_file()):
        sys.exit(f"engine binaries not found in {REPO}: run `make windows` "
                 "(Windows) or rerun with --build-engine")

    run([find_pyinstaller(), "--noconfirm", "--distpath", "gui/dist",
         "--workpath", "gui/build", "gui/studio.spec"])

    bundle = REPO / "dist" / f"fsmp-kmc-{args.name}-{platform_tag()}"
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
    for folder in ("configs", "models", "cells"):
        shutil.copytree(REPO / folder, bundle / folder)
    (bundle / "forcefields").mkdir()
    shutil.copy2(REPO / "forcefields" / "readme.txt", bundle / "forcefields")
    shutil.copy2(REPO / ".github" / "release_readme.txt",
                 bundle / "README.txt")

    kind = "zip" if os.name == "nt" else "gztar"
    archive = shutil.make_archive(str(bundle), kind,
                                  root_dir=bundle.parent,
                                  base_dir=bundle.name)
    print("bundle: ", bundle)
    print("archive:", archive)


if __name__ == "__main__":
    main()
