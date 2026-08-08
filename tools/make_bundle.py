#!/usr/bin/env python3
"""Assemble the release archives: the full bundle (the Studio frozen with
PyInstaller, the engine binaries and the bundled data) and the -cli archive
(the same minus the Studio, a few megabytes for command-line use), laid out
exactly like a published release.

    python3 tools/make_bundle.py [--name v0.5.0] [--build-engine] [--cli]

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
from datetime import datetime, timezone
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


# Everything an archive must carry whatever else is in it: the engine, the
# papers that let it be redistributed, and the data the demonstration runs on
COMMON_REQUIRED = [
    "fsmp" + EXE, "pack" + EXE,
    "README.txt", "LICENSE", "THIRD_PARTY_NOTICES.md", "BUILD_INFO.txt",
    "forcefields/readme.txt",
    "configs/tma_quickstart_demo.txt",
    "samples/models/trimesic_acid.xyz",
    "samples/projects/TMA_quickstart/project.json",
]


def verify(bundle: Path, required: list) -> None:
    """Refuse to package an archive that lost something a user needs."""
    missing = [name for name in required if not (bundle / name).exists()]
    # the demonstration is useless without its grid; the name is free to
    # change, the folder may not be empty
    if not list((bundle / "samples" / "potentials").glob("*.v2.bin")):
        missing.append("samples/potentials/*.v2.bin")
    if missing:
        sys.exit(f"{bundle.name} is incomplete, refusing to package:\n  "
                 + "\n  ".join(missing))
    print(f"verified {len(required) + 1} required entries in {bundle.name}")


def git_commit() -> str:
    """The commit an archive was built from, marked when the tree carried
    uncommitted changes. Unknown outside a git checkout."""
    try:
        sha = subprocess.run(["git", "rev-parse", "--short=12", "HEAD"],
                             cwd=REPO, capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO,
                               capture_output=True, text=True,
                               check=True).stdout.strip()
        return sha + (" (with uncommitted changes)" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown (not a git checkout)"


def engine_version(engine: Path) -> str:
    try:
        return subprocess.run([str(engine), "--version"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as e:
        return f"could not be asked ({e})"


def write_build_info(bundle: Path, name: str, engine: Path,
                     studio: bool) -> None:
    """A short record of where this archive came from, so a bug report can
    name the build instead of guessing at it."""
    lines = [
        f"FSMP-kMC {project_version()}",
        f"archive:   {bundle.name}",
        f"contents:  {'Studio, engine and data' if studio else 'engine and data, no Studio'}",
        f"built:     {datetime.now(timezone.utc):%Y-%m-%d %H:%M:%S} UTC",
        f"platform:  {platform_tag()} ({platform.platform()})",
        f"commit:    {git_commit()}",
        f"engine:    {engine_version(engine)}",
    ]
    if studio:
        lines.append(f"python:    {platform.python_version()}")
        frozen = []
        for package in ("PySide6", "numpy", "rdkit", "pillow", "pyinstaller"):
            try:
                frozen.append(f"{package} {metadata.version(package)}")
            except metadata.PackageNotFoundError:
                pass
        lines.append("frozen:    " + ", ".join(frozen))
    (bundle / "BUILD_INFO.txt").write_text("\n".join(lines) + "\n",
                                           encoding="utf-8")


def copy_shared_data(bundle: Path, engine: Path, pack: Path) -> None:
    """The parts both archives carry: the binaries, the examples and the
    empty forcefields folder the user drops downloads into."""
    shutil.copy2(engine, bundle)
    shutil.copy2(pack, bundle)
    # configs are the engine's example parameter files; samples carries the
    # example molecule models, unit cells, the demonstration potential and
    # the ready-to-open Studio projects
    for folder in ("configs", "samples"):
        shutil.copytree(REPO / folder, bundle / folder)
    (bundle / "forcefields").mkdir()
    shutil.copy2(REPO / "forcefields" / "readme.txt", bundle / "forcefields")


def build_cli_archive(name: str, engine: Path, pack: Path) -> str:
    """The engine on its own: a few megabytes for someone who runs from the
    command line and has no use for a hundred megabytes of Qt and Python.
    The staging folder is removed once packed, so dist/ keeps exactly one
    unpacked bundle (the release workflow self-tests it by globbing)."""
    stage = REPO / "dist" / f"fsmp-kmc-{name}-{platform_tag()}-cli"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)

    copy_shared_data(stage, engine, pack)
    shutil.copy2(REPO / ".github" / "release_readme_cli.txt",
                 stage / "README.txt")
    shutil.copy2(REPO / "LICENSE", stage / "LICENSE")
    shutil.copy2(REPO / "THIRD_PARTY_NOTICES.md", stage)
    write_build_info(stage, name, engine, studio=False)
    verify(stage, COMMON_REQUIRED)

    kind = "zip" if os.name == "nt" else "gztar"
    archive = shutil.make_archive(str(stage), kind, root_dir=stage.parent,
                                  base_dir=stage.name)
    shutil.rmtree(stage)
    return archive


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default=None,
                        help="version part of the bundle name (e.g. v0.5.0);"
                             " defaults to version.h")
    parser.add_argument("--build-engine", action="store_true",
                        help="compile the engine here with the release flags")
    parser.add_argument("--cli", action="store_true",
                        help="only the command-line archive (no PyInstaller)")
    args = parser.parse_args()
    name = args.name or ("v" + project_version())

    if args.build_engine:
        build_engine()
    engine = REPO / ("fsmp" + EXE)
    pack = REPO / ("pack" + EXE)
    if not (engine.is_file() and pack.is_file()):
        sys.exit(f"engine binaries not found in {REPO}: run `make windows` "
                 "(Windows) or rerun with --build-engine")

    if args.cli:
        print("archive:", build_cli_archive(name, engine, pack))
        return

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
    copy_shared_data(bundle, engine, pack)
    shutil.copy2(REPO / ".github" / "release_readme.txt",
                 bundle / "README.txt")
    copy_licenses(bundle)
    write_build_info(bundle, name, engine, studio=True)
    verify(bundle, COMMON_REQUIRED + [studio_name(), "licenses/LGPL-3.0.txt"])

    kind = "zip" if os.name == "nt" else "gztar"
    archive = shutil.make_archive(str(bundle), kind,
                                  root_dir=bundle.parent,
                                  base_dir=bundle.name)
    print("bundle: ", bundle)
    print("archive:", archive)
    print("archive:", build_cli_archive(name, engine, pack))


if __name__ == "__main__":
    main()
