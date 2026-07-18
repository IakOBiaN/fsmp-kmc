"""Molecule model and xyz file input/output.

The geometry is deliberately plain: element and coordinates only.
Charges and other interaction parameters belong to the potential model,
not to the molecule file.
"""

import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .elements import BOND_TOLERANCE, covalent_radius, normalize_symbol


@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float = 0.0


def bonded_pairs(atoms: list["Atom"]) -> list[tuple[int, int]]:
    """Index pairs close enough to be bonded, by the covalent-radius rule.
    This is the single source of truth for connectivity: the editor draws these
    exact bonds, and the optimizer and potential generator perceive from them,
    so what the user sees on the canvas is what gets built."""
    pairs = []
    for i in range(len(atoms)):
        ai = atoms[i]
        ri = covalent_radius(ai.element)
        for j in range(i + 1, len(atoms)):
            aj = atoms[j]
            d = math.dist((ai.x, ai.y, ai.z), (aj.x, aj.y, aj.z))
            if 1e-6 < d < BOND_TOLERANCE * (ri + covalent_radius(aj.element)):
                pairs.append((i, j))
    return pairs


def centroid(atoms: list["Atom"]) -> tuple[float, float]:
    """The unweighted centre of the atom positions in the plane. The engine
    rotates a molecule about the origin (molecule_model.h), so a finished
    model normally has this point at (0, 0)."""
    n = len(atoms)
    return (sum(a.x for a in atoms) / n, sum(a.y for a in atoms) / n)


def translated(atoms: list["Atom"], dx: float, dy: float) -> list["Atom"]:
    """A copy of the geometry shifted in the plane (z is untouched)."""
    return [Atom(a.element, a.x + dx, a.y + dy, a.z) for a in atoms]


def rotated(atoms: list["Atom"], angle_deg: float) -> list["Atom"]:
    """A copy of the geometry rotated rigidly about the origin,
    counterclockwise for a positive angle: the same sense as the orientation
    angle of a molecule placed in a cell."""
    rad = math.radians(angle_deg)
    c, s = math.cos(rad), math.sin(rad)
    return [Atom(a.element, a.x * c - a.y * s, a.x * s + a.y * c, a.z)
            for a in atoms]


def aimed_at_x(atoms: list["Atom"], index: int) -> list["Atom"]:
    """A copy of the geometry rotated rigidly about the origin so that atom
    `index` lands on the positive x axis, which is the reference direction of
    the zero orientation angle. Raises ValueError for an atom at the origin,
    which defines no direction."""
    a = atoms[index]
    r = math.hypot(a.x, a.y)
    if r < 1e-6:
        raise ValueError("the chosen atom sits at the rotation centre; "
                         "pick one away from the origin")
    out = rotated(atoms, -math.degrees(math.atan2(a.y, a.x)))
    out[index] = Atom(a.element, r, 0.0, a.z)   # exactly on the axis
    return out


def connected_components(atoms: list["Atom"]) -> list[list[int]]:
    """Groups of atom indices joined through bonds (see bonded_pairs). A whole
    molecule is a single group; more than one means the geometry has detached
    atoms or parts that are not bonded to the rest."""
    parent = list(range(len(atoms)))

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:       # path compression
            parent[x], x = root, parent[x]
        return root

    for i, j in bonded_pairs(atoms):
        parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in range(len(atoms)):
        groups.setdefault(find(i), []).append(i)
    return list(groups.values())


class Molecule:
    def __init__(self, atoms: list[Atom] | None = None, comment: str = ""):
        self.atoms: list[Atom] = atoms if atoms is not None else []
        self.comment = comment

    def formula(self) -> str:
        """Hill order: C, H, then the rest alphabetically."""
        counts = Counter(a.element for a in self.atoms)
        parts = []
        for el in ["C", "H"] + sorted(k for k in counts if k not in ("C", "H")):
            n = counts.get(el, 0)
            if n:
                parts.append(el + (str(n) if n > 1 else ""))
        return "".join(parts)

    @staticmethod
    def load_xyz(path: str | Path) -> "Molecule":
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        if not lines:
            raise ValueError("empty file")
        try:
            count = int(lines[0].split()[0])
        except (ValueError, IndexError):
            raise ValueError("line 1: expected the number of atoms")
        comment = lines[1] if len(lines) > 1 else ""
        atoms = []
        for i in range(count):
            ln = 2 + i
            if ln >= len(lines):
                raise ValueError(f"expected {count} atoms, file ends after {i}")
            cols = lines[ln].split()
            if len(cols) < 4:
                raise ValueError(f"line {ln + 1}: expected 'element x y z'")
            el = normalize_symbol(cols[0])
            try:
                x, y, z = float(cols[1]), float(cols[2]), float(cols[3])
            except ValueError:
                raise ValueError(f"line {ln + 1}: bad coordinate")
            atoms.append(Atom(el, x, y, z))  # extra columns are ignored
        return Molecule(atoms, comment)

    def save_xyz(self, path: str | Path) -> None:
        lines = [str(len(self.atoms))]
        lines.append(self.comment if self.comment else self.formula())
        for a in self.atoms:
            lines.append(f"{a.element:<2} {a.x:15.8f} {a.y:15.8f} {a.z:15.8f}")
        Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
