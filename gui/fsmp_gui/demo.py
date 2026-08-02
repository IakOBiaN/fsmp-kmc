"""The one-click demonstration: a copy of the bundled TMA quickstart project
that a first-time user can run immediately.

The shipped copy in samples/projects/ is reference data. Running it in place
would write run folders into it, and an archive unpacked somewhere read-only
could not be run at all, so the demonstration is installed into the user's
documents folder first. The copy carries its own potential, which keeps it
self-contained: it survives being moved, and it holds no path that only
exists on the machine that made it.

No Qt in this module: everything is testable standalone.
"""

import shutil
from pathlib import Path

from .engine import app_root
from .project import Project

PROJECT_NAME = "TMA_quickstart"
GRID_NAME = "TMA_simple_2020_coarse_demo.v2.bin"
INSTALL_FOLDER = "FSMP-kMC"
POTENTIAL_DIR = "potentials"


class DemoError(Exception):
    pass


def source_project() -> Path:
    return app_root() / "samples" / "projects" / PROJECT_NAME


def source_grid() -> Path:
    return app_root() / "samples" / "potentials" / GRID_NAME


def available() -> bool:
    """True when this build actually carries the demonstration data."""
    return (source_project() / "project.json").is_file() and source_grid().is_file()


def install(parent: Path) -> tuple[Project, bool]:
    """Install the demonstration into parent/TMA_quickstart and open it.
    Returns (project, created); created is False when the copy was already
    there, so a second click reopens the same project with its earlier runs
    instead of duplicating it."""
    target = Path(parent) / PROJECT_NAME
    if (target / "project.json").is_file():
        return Project.open(target), False
    if not available():
        raise DemoError("this build does not carry the demonstration project "
                        f"(expected {source_project()})")
    if target.exists() and any(target.iterdir()):
        raise DemoError(f"{target} exists and is not a project folder")

    try:
        shutil.copytree(source_project(), target, dirs_exist_ok=True)
        grid = target / POTENTIAL_DIR / GRID_NAME
        grid.parent.mkdir(exist_ok=True)
        shutil.copy2(source_grid(), grid)
    except OSError as e:
        raise DemoError(f"cannot install the demonstration: {e}")

    # the shipped manifest points at the shared samples/potentials folder;
    # the copy owns its grid, so the path becomes project-relative
    project = Project.open(target)
    project.set_potential(Path(GRID_NAME).name.split(".v2")[0], grid)
    return project, True
