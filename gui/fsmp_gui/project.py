"""Project model: a folder with a project.json manifest and data subfolders.

The manifest keeps lightweight references (names and relative paths); the
data itself lives in files inside the project folder, so a project stays
usable without the GUI.
"""

import json
import re
from datetime import datetime, timezone
from pathlib import Path

MANIFEST = "project.json"
FORMAT = 1


class ProjectError(Exception):
    pass


def safe_filename(name: str) -> str:
    s = re.sub(r'[<>:"/\\|?*\s]+', "_", name.strip()).strip("_.")
    return s or "unnamed"


class Project:
    def __init__(self, root: Path, manifest: dict):
        self.root = Path(root)
        self.manifest = manifest

    @property
    def name(self) -> str:
        return self.manifest.get("name", self.root.name)

    @property
    def molecule(self) -> dict | None:
        """The project molecule ({"name", "file"}) or None. A project holds
        exactly one molecule."""
        return self.manifest.get("molecule")

    @classmethod
    def create(cls, root: str | Path, name: str) -> "Project":
        root = Path(root)
        if root.exists() and any(root.iterdir()):
            raise ProjectError(f"folder is not empty: {root}")
        (root / "molecules").mkdir(parents=True, exist_ok=True)
        manifest = {
            "format": FORMAT,
            "name": name,
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "molecule": None,
            "potentials": [],
            "unit_cell": None,
            "simulation_cell": None,
            "simulation": None,
        }
        project = cls(root, manifest)
        project.save()
        return project

    @classmethod
    def open(cls, root: str | Path) -> "Project":
        root = Path(root)
        path = root / MANIFEST
        if not path.is_file():
            raise ProjectError(f"not a project folder (no {MANIFEST}): {root}")
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as e:
            raise ProjectError(f"broken {MANIFEST}: {e}")
        if manifest.get("format") != FORMAT:
            raise ProjectError(f"unsupported project format: {manifest.get('format')!r}")
        return cls(root, manifest)

    def save(self) -> None:
        path = self.root / MANIFEST
        path.write_text(json.dumps(self.manifest, indent=2) + "\n", encoding="utf-8")

    # -- the project molecule ----------------------------------------------

    def molecule_path(self, entry: dict | None = None) -> Path:
        entry = entry if entry is not None else self.molecule
        return self.root / entry["file"]

    def set_molecule(self, name: str, molecule) -> dict:
        """Save the molecule and make it the project molecule, replacing
        the previous one (and its file) if any."""
        (self.root / "molecules").mkdir(exist_ok=True)
        rel = f"molecules/{safe_filename(name)}.xyz"
        molecule.save_xyz(self.root / rel)
        old = self.molecule
        if old is not None and old["file"] != rel:
            old_path = self.root / old["file"]
            if old_path.is_file():
                old_path.unlink()
        self.manifest["molecule"] = {"name": name, "file": rel}
        self.save()
        return self.manifest["molecule"]

    def clear_molecule(self) -> None:
        entry = self.molecule
        if entry is None:
            return
        path = self.molecule_path(entry)
        if path.is_file():
            path.unlink()
        self.manifest["molecule"] = None
        self.save()
