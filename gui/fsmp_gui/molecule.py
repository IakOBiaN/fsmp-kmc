"""Molecule model and xyz file input/output.

The geometry is deliberately plain: element and coordinates only.
Charges and other interaction parameters belong to the potential model,
not to the molecule file.
"""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from .elements import normalize_symbol


@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float = 0.0


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
