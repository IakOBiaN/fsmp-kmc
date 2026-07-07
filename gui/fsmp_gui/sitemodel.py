"""Coarse-grained "site" model of a molecule: a set of interaction sites,
each a Lennard-Jones center and/or a point charge.

Stored in an xyz-like text format (.site), one site per line:

    <n_sites>
    <comment / model name>
    LABEL   x   y   z   q   epsilon   sigma

x, y, z in angstroms; q in elementary charges; epsilon in kelvin (eps/k_B);
sigma in angstroms. A pure charge has epsilon = 0, a pure LJ center q = 0.
"""

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

SUFFIX = ".site"
_HEADER = "LABEL x y z q(e) epsilon(K) sigma(A)"


@dataclass
class Site:
    label: str
    x: float
    y: float
    z: float = 0.0
    q: float = 0.0
    epsilon: float = 0.0
    sigma: float = 0.0

    @property
    def is_lj(self) -> bool:
        return self.epsilon > 0.0 and self.sigma > 0.0

    @property
    def is_charge(self) -> bool:
        return self.q != 0.0


class SiteModel:
    def __init__(self, sites: list[Site] | None = None, comment: str = ""):
        self.sites: list[Site] = sites if sites is not None else []
        self.comment = comment

    def total_charge(self) -> float:
        return sum(s.q for s in self.sites)

    def summary(self) -> str:
        n_lj = sum(1 for s in self.sites if s.is_lj)
        n_q = sum(1 for s in self.sites if s.is_charge)
        return f"{len(self.sites)} sites ({n_lj} LJ, {n_q} charged)"

    def label_counts(self) -> dict[str, int]:
        return dict(Counter(s.label for s in self.sites))

    @staticmethod
    def load(path: str | Path) -> "SiteModel":
        lines = Path(path).read_text(encoding="utf-8").splitlines()
        if not lines:
            raise ValueError("empty file")
        try:
            count = int(lines[0].split()[0])
        except (ValueError, IndexError):
            raise ValueError("line 1: expected the number of sites")
        comment = lines[1] if len(lines) > 1 else ""
        sites = []
        for i in range(count):
            ln = 2 + i
            if ln >= len(lines):
                raise ValueError(f"expected {count} sites, file ends after {i}")
            cols = lines[ln].split()
            if len(cols) < 7:
                raise ValueError(f"line {ln + 1}: expected "
                                 "'label x y z q epsilon sigma'")
            try:
                x, y, z, q, eps, sig = (float(c) for c in cols[1:7])
            except ValueError:
                raise ValueError(f"line {ln + 1}: bad number")
            sites.append(Site(cols[0], x, y, z, q, eps, sig))
        return SiteModel(sites, comment)

    def save(self, path: str | Path) -> None:
        lines = [str(len(self.sites)), self.comment or _HEADER]
        for s in self.sites:
            lines.append(f"{s.label:<6} {s.x:13.6f} {s.y:13.6f} {s.z:13.6f} "
                         f"{s.q:10.5f} {s.epsilon:12.4f} {s.sigma:10.5f}")
        Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
