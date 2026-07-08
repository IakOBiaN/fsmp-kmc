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
    def model(self) -> dict | None:
        """The project molecule model, or None. A project holds exactly one.
        {"kind": "atomistic"|"site", "name", "file"}."""
        return self.manifest.get("model")

    @property
    def model_kind(self) -> str | None:
        m = self.model
        return m["kind"] if m else None

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
            "model": None,
            "potential": None,
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

    # -- the project molecule model ----------------------------------------

    def model_path(self, entry: dict | None = None) -> Path:
        entry = entry if entry is not None else self.model
        return self.root / entry["file"]

    def _set_model(self, kind: str, name: str, suffix: str, save_fn) -> dict:
        """Replace the project model (of either kind), writing its file and
        deleting the previous one."""
        (self.root / "molecules").mkdir(exist_ok=True)
        rel = f"molecules/{safe_filename(name)}{suffix}"
        save_fn(self.root / rel)
        old = self.model
        if old is not None and old["file"] != rel:
            old_path = self.root / old["file"]
            if old_path.is_file():
                old_path.unlink()
        self.manifest["model"] = {"kind": kind, "name": name, "file": rel}
        self.save()
        return self.manifest["model"]

    def set_atomistic(self, name: str, molecule) -> dict:
        return self._set_model("atomistic", name, ".xyz", molecule.save_xyz)

    def set_site_model(self, name: str, site_model) -> dict:
        return self._set_model("site", name, ".site", site_model.save)

    def clear_model(self) -> None:
        entry = self.model
        if entry is None:
            return
        path = self.model_path(entry)
        if path.is_file():
            path.unlink()
        self.manifest["model"] = None
        self.save()

    # -- the project potential ----------------------------------------------

    @property
    def potential(self) -> dict | None:
        """{"name", "path"} or None. The file is referenced, not copied:
        forcefields are often hundreds of megabytes."""
        return self.manifest.get("potential")

    def potential_path(self) -> Path:
        entry = self.potential
        p = Path(entry["path"])
        return p if p.is_absolute() else self.root / p

    def set_potential(self, name: str, path: str | Path) -> dict:
        p = Path(path).resolve()
        try:
            stored = p.relative_to(self.root.resolve()).as_posix()
        except ValueError:
            stored = str(p)
        self.manifest["potential"] = {"name": name, "path": stored}
        self.save()
        return self.manifest["potential"]

    def clear_potential(self) -> None:
        self.manifest["potential"] = None
        self.save()
