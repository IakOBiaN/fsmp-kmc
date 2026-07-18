"""Generate a pair potential from a molecule model.

Two rigid copies of the molecule are swept over a grid of (r, alpha1, alpha2):
molecule A sits at the origin with orientation alpha1, molecule B at (r, 0)
with orientation alpha2, both rotating about the origin (the molecule centre,
per molecule_model.h). An energy backend turns each (r, angle grid) into a raw
slab in J/mol; the slab is capped, exchange-symmetrised and folded exactly as
tools/pack_forcefield.cpp does, then written in the binary v2 format.

Two backends share this machinery:

- SiteBackend: Lennard-Jones + Coulomb between the coarse-grained sites of a
  site model (the original, released generator);
- MMFFBackend: the full atomistic molecule scored with MMFF94 (buffered-14-7
  van der Waals + buffered Coulomb) on RDKit's typed parameters, so a project
  with only an atomistic model can still make a potential.
"""

import math
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .sitemodel import SiteModel

R = 8.314462618                     # J/(mol K)
Ke = 8.9875517923e19               # J*A/C^2
Na = 6.02214076e23
E_CHARGE = 1.602176634e-19
CAP_JMOL = 1.0e4 * 4184.0          # 1e4 kcal/mol, matches the reference grid
KCAL_TO_JMOL = 4184.0
MMFF_ELE = 332.0716                # kcal*A/(mol*e^2), MMFF electrostatic constant
MMFF_DELTA = 0.05                  # A, MMFF electrostatic buffering

HEADER_BYTES = 64
MAGIC = b"FSMP"
VERSION = 2


@dataclass
class GridSpec:
    r_min: float
    r_max: float
    dr: float
    da: float          # angular step, degrees
    fold_deg: float    # 360 = no folding
    use_float: bool = True

    @property
    def n_dist(self) -> int:
        return int(round((self.r_max - self.r_min) / self.dr)) + 1


class GenerationError(Exception):
    pass


def _rotated(coords: np.ndarray, angles_rad: np.ndarray) -> np.ndarray:
    """coords (n,2) local -> (len(angles), n, 2) rotated about the origin."""
    c, s = np.cos(angles_rad), np.sin(angles_rad)
    x, y = coords[:, 0], coords[:, 1]
    rx = np.outer(c, x) - np.outer(s, y)
    ry = np.outer(s, x) + np.outer(c, y)
    return np.stack([rx, ry], axis=-1)


def _prepare(model: SiteModel):
    lj = [(s.x, s.y, s.sigma, s.epsilon * R, s.r0) for s in model.sites if s.is_lj]
    ch = [(s.x, s.y, s.q) for s in model.sites if s.q != 0.0]
    if not lj and not ch:
        raise GenerationError("the site model has no LJ centers and no charges")
    return lj, ch


def _slab(r: float, ang: np.ndarray, lj, ch) -> np.ndarray:
    """Raw energy slab (na x na) at centre-centre distance r, J/mol.
    Rows index alpha1 (molecule A at origin), columns alpha2 (molecule B at (r,0))."""
    na = len(ang)
    e = np.zeros((na, na))
    # LJ: combine pairs with Lorentz-Berthelot + averaged hard core
    if lj:
        A = _rotated(np.array([(x, y) for x, y, *_ in lj]), ang)          # (na, nl, 2)
        B = _rotated(np.array([(x, y) for x, y, *_ in lj]), ang)
        B = B + np.array([r, 0.0])
        for i, (_, _, si, ei, r0i) in enumerate(lj):
            for j, (_, _, sj, ej, r0j) in enumerate(lj):
                sig = 0.5 * (si + sj)
                eps = math.sqrt(ei * ej)
                r0 = 0.5 * (r0i + r0j)
                dx = A[:, i, 0][:, None] - B[:, j, 0][None, :]
                dy = A[:, i, 1][:, None] - B[:, j, 1][None, :]
                d = np.sqrt(dx * dx + dy * dy)
                x6 = ((sig - r0) / np.maximum(d - r0, 1e-9)) ** 6
                e += 4.0 * eps * (x6 * (x6 - 1.0))
    # Coulomb
    if ch:
        qa = np.array([q for *_, q in ch])
        A = _rotated(np.array([(x, y) for x, y, _ in ch]), ang)
        B = _rotated(np.array([(x, y) for x, y, _ in ch]), ang) + np.array([r, 0.0])
        const = Ke * Na * (E_CHARGE ** 2)
        for i, qi in enumerate(qa):
            for j, qj in enumerate(qa):
                dx = A[:, i, 0][:, None] - B[:, j, 0][None, :]
                dy = A[:, i, 1][:, None] - B[:, j, 1][None, :]
                d = np.sqrt(dx * dx + dy * dy)
                e += const * qi * qj / np.maximum(d, 1e-9)
    return e


def _fold(slab: np.ndarray, spec: GridSpec):
    """Exchange-symmetrise (max of the two partners) then fold-average,
    mirroring tools/pack_forcefield.cpp. Returns the output slab."""
    na = slab.shape[0]                        # raw angles 0..360 inclusive
    full = na - 1                             # 360 in index units
    half = full // 2
    # exchange symmetry: U(t1,t2) = U(t2+180, t1+180), keep the max.
    # partner[j,k] = slab[(k+half)%full, (j+half)%full]
    idx = (np.arange(na) + half) % full
    partner = slab[np.ix_(idx, idx)].T
    out = np.maximum(slab, partner)
    if spec.fold_deg >= 359.999:
        return out
    fold_p = int(round(spec.fold_deg / spec.da))
    copies = int(round(360.0 / spec.fold_deg))
    on = fold_p + 1
    folded = np.zeros((on, on))
    for a in range(copies):
        for b in range(copies):
            folded += out[a * fold_p: a * fold_p + on, b * fold_p: b * fold_p + on]
    return folded / (copies * copies)


class SiteBackend:
    """Lennard-Jones + Coulomb energy between the sites of a site model."""

    def __init__(self, model: SiteModel):
        self.lj, self.ch = _prepare(model)

    def slab(self, r: float, ang: np.ndarray) -> np.ndarray:
        return _slab(r, ang, self.lj, self.ch)


class MMFFBackend:
    """MMFF94 energy between two rigid copies of an atomistic molecule.

    RDKit's own force field does not score the interaction between two disjoint
    fragments, so the intermolecular energy is summed here directly from the
    typed parameters: buffered-14-7 van der Waals plus a buffered Coulomb term
    on the MMFF partial charges. Result is J/mol, matching SiteBackend."""

    def __init__(self, molecule, net_charge: int = 0):
        from . import mmff
        p = mmff.mmff_pair_params(molecule, net_charge)
        self.xy = np.ascontiguousarray(p.xy, dtype=float)   # (n, 2)
        self.q = np.ascontiguousarray(p.q, dtype=float)     # (n,)
        self.Rstar = np.ascontiguousarray(p.Rstar, dtype=float)
        self.eps = np.ascontiguousarray(p.eps, dtype=float)
        self.R7 = self.Rstar ** 7

    def slab(self, r: float, ang: np.ndarray) -> np.ndarray:
        na = len(ang)
        n = self.xy.shape[0]
        A = _rotated(self.xy, ang)                       # (na, n, 2)
        B = _rotated(self.xy, ang) + np.array([r, 0.0])
        e = np.zeros((na, na))
        for i in range(n):
            axi = A[:, i, 0][:, None]
            ayi = A[:, i, 1][:, None]
            qi = self.q[i]
            for j in range(n):
                dx = axi - B[:, j, 0][None, :]
                dy = ayi - B[:, j, 1][None, :]
                d = np.sqrt(dx * dx + dy * dy)
                np.maximum(d, 1e-6, out=d)
                Rs, ep, R7 = self.Rstar[i, j], self.eps[i, j], self.R7[i, j]
                t = (1.07 * Rs) / (d + 0.07 * Rs)
                evdw = ep * t ** 7 * (1.12 * R7 / (d ** 7 + 0.12 * R7) - 2.0)
                e += evdw + MMFF_ELE * qi * self.q[j] / (d + MMFF_DELTA)
        return e * KCAL_TO_JMOL


def generate(backend, spec: GridSpec, out_path: str | Path,
             progress=None, cancel=None) -> bool:
    """Write the v2 potential. `backend.slab(r, ang)` returns the raw J/mol
    energy matrix; a bare SiteModel is accepted too, for backward compatibility.
    `progress(done, total)` is called per distance row. `cancel()` returning
    True aborts, removes the partial file and makes this return False; a fully
    written grid returns True."""
    if isinstance(backend, SiteModel):
        backend = SiteBackend(backend)
    na_raw = int(round(360.0 / spec.da)) + 1
    ang = np.deg2rad(np.arange(na_raw) * spec.da)
    if spec.fold_deg < 359.999:
        fold_p = int(round(spec.fold_deg / spec.da))
        if abs(fold_p * spec.da - spec.fold_deg) > 1e-9 or 360 % round(spec.fold_deg):
            raise GenerationError("fold angle must divide 360 and be a multiple of da")
        out_n = fold_p + 1
        out_fold = spec.fold_deg
    else:
        out_n = na_raw
        out_fold = 360.0

    n_dist = spec.n_dist
    dtype = np.float32 if spec.use_float else np.float64
    out_path = Path(out_path)
    with out_path.open("wb") as f:
        header = (MAGIC
                  + struct.pack("<4I", VERSION, 1 if spec.use_float else 0,
                                n_dist, out_n)
                  + struct.pack("<4d", spec.r_min, spec.dr, spec.da, out_fold))
        f.write(header + b"\x00" * (HEADER_BYTES - len(header)))
        for i in range(n_dist):
            if cancel is not None and cancel():
                f.close()
                out_path.unlink(missing_ok=True)
                return False
            r = spec.r_min + i * spec.dr
            slab = backend.slab(r, ang)
            np.clip(slab, -CAP_JMOL, CAP_JMOL, out=slab)
            folded = _fold(slab, spec)
            folded.astype(dtype).tofile(f)
            if progress is not None:
                progress(i + 1, n_dist)
    return True
