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
    def atomistic(self) -> dict | None:
        """The atomistic model {"name", "file"} or None. Used for visualization
        (and, on its own, for future atomistic potential generation)."""
        return self.manifest.get("atomistic")

    @property
    def site(self) -> dict | None:
        """The site (charge) model {"name", "file"} or None. Used to generate
        the pair potential."""
        return self.manifest.get("site")

    @property
    def generation_kind(self) -> str | None:
        """Which model drives potential generation: the site model if present,
        otherwise the atomistic one."""
        if self.site is not None:
            return "site"
        if self.atomistic is not None:
            return "atomistic"
        return None

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
            "atomistic": None,
            "site": None,
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

    # -- the project molecule models ---------------------------------------
    # A project holds up to two independent models: an atomistic one (for
    # visualization) and a site one (for the pair potential).

    def model_path(self, entry: dict) -> Path:
        return self.root / entry["file"]

    def _save_model(self, key: str, name: str, suffix: str, save_fn) -> dict:
        (self.root / "molecules").mkdir(exist_ok=True)
        rel = f"molecules/{safe_filename(name)}{suffix}"
        save_fn(self.root / rel)
        old = self.manifest.get(key)
        if old is not None and old["file"] != rel:
            old_path = self.root / old["file"]
            if old_path.is_file():
                old_path.unlink()
        self.manifest[key] = {"name": name, "file": rel}
        self.save()
        return self.manifest[key]

    def _clear_model(self, key: str) -> None:
        old = self.manifest.get(key)
        if old is not None:
            path = self.root / old["file"]
            if path.is_file():
                path.unlink()
        self.manifest[key] = None
        self.save()

    def set_atomistic(self, name: str, molecule) -> dict:
        return self._save_model("atomistic", name, ".xyz", molecule.save_xyz)

    def set_site(self, name: str, site_model) -> dict:
        return self._save_model("site", name, ".site", site_model.save)

    def clear_atomistic(self) -> None:
        self._clear_model("atomistic")

    def clear_site(self) -> None:
        self._clear_model("site")

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

    # -- the rough unit cell -------------------------------------------------

    @property
    def unit_cell(self) -> dict | None:
        """{"cell_x", "cell_y", "molecules": [{"x","y","phi"}, ...]} or None."""
        return self.manifest.get("unit_cell")

    def set_unit_cell(self, cell_x: float, cell_y: float,
                      molecules: list[dict]) -> dict:
        self.manifest["unit_cell"] = {"cell_x": cell_x, "cell_y": cell_y,
                                      "molecules": list(molecules)}
        self.save()
        return self.manifest["unit_cell"]

    def clear_unit_cell(self) -> None:
        self.manifest["unit_cell"] = None
        self.save()

    # -- the elongated simulation cell ----------------------------------------

    @property
    def simulation_cell(self) -> dict | None:
        """Settings of the elongated cell (uc_in_x, uc_in_y, free_space,
        lambdam, temperature_in_transition_zone, preview values) or None."""
        return self.manifest.get("simulation_cell")

    def set_simulation_cell(self, settings: dict) -> dict:
        self.manifest["simulation_cell"] = dict(settings)
        self.save()
        return self.manifest["simulation_cell"]

    # -- the run form defaults -------------------------------------------------

    @property
    def simulation(self) -> dict | None:
        """The last used Run-tab form values, so a new run starts from the
        previous settings."""
        return self.manifest.get("simulation")

    def set_simulation(self, form: dict) -> dict:
        self.manifest["simulation"] = dict(form)
        self.save()
        return self.manifest["simulation"]
