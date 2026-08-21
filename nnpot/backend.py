import numpy as np

from .geometry import (ATOMIC_NUMBER, dimer_coords, orbit_map, pair_distances,
                       rotational_order)

CAP_JMOL = 1.0e4 * 4184.0
KE = 8.9875517923e19 * 6.02214076e23 * (1.602176634e-19) ** 2
ATOM_BUDGET = 90000


def _switch(d, lo, hi):
    t = np.clip((d - lo) / max(hi - lo, 1e-12), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


class NNBackend:

    def __init__(self, elements, xy, calculator, da, charges=None, order=None,
                 wall=1.3, cap=CAP_JMOL, lr_window=1.0, chunk=None,
                 progress=None):
        self.elements = list(elements)
        self.xy = np.ascontiguousarray(xy, dtype=float)
        self.calc = calculator
        self.da = float(da)
        self.wall = float(wall)
        self.cap = float(cap)
        self.lr_window = float(lr_window)
        self.chunk = int(chunk) if chunk else max(1, ATOM_BUDGET // (2 * len(self.xy)))
        self.progress = progress
        self.charges = None if charges is None else np.asarray(charges, float)
        self.order = (rotational_order(self.xy, self.elements)
                      if order is None else int(order))
        self.radius = float(np.hypot(self.xy[:, 0], self.xy[:, 1]).max())
        missing = sorted({e for e in self.elements if e not in ATOMIC_NUMBER})
        if missing:
            raise ValueError("no atomic number for " + ", ".join(missing))
        self.z = np.array([ATOMIC_NUMBER[e] for e in self.elements], dtype=int)
        self.z2 = np.concatenate([self.z, self.z])
        self.cutoff = getattr(calculator, "cutoff", None)
        mono = np.concatenate([self.xy, np.zeros((len(self.xy), 1))], axis=1)
        self.e_mono = float(np.asarray(calculator.energies(self.z, mono[None]))[0])
        self.stats = {"evaluated": 0, "walled": 0, "rows": 0, "rows_beyond_reach": 0}
        self._orbits = {}

    def _orbit(self, p):
        if p not in self._orbits:
            self._orbits[p] = orbit_map(p, self.da)
        return self._orbits[p]

    def _long_range(self, r, t1, t2):
        d = pair_distances(self.xy, r, t1, t2)
        qq = self.charges[:, None] * self.charges[None, :]
        w = _switch(d, self.cutoff - self.lr_window, self.cutoff)
        return (KE * qq * w / np.maximum(d, 1e-9)).sum(axis=(1, 2))

    def _energies(self, r, t1, t2):
        total = len(t1)
        out = np.zeros(total)
        beyond = self.cutoff is not None and r - 2.0 * self.radius >= self.cutoff
        if beyond:
            self.stats["rows_beyond_reach"] += 1
        for start in range(0, total, self.chunk):
            stop = min(start + self.chunk, total)
            a1, a2 = t1[start:stop], t2[start:stop]
            d = pair_distances(self.xy, r, a1, a2)
            dmin = d.reshape(len(a1), -1).min(axis=1)
            walled = dmin < self.wall
            chunk_out = np.zeros(len(a1))
            chunk_out[walled] = self.cap
            live = ~walled
            if not beyond and live.any():
                coords = dimer_coords(self.xy, r, a1[live], a2[live])
                energy = np.asarray(self.calc.energies(self.z2, coords), dtype=float)
                chunk_out[live] = energy - 2.0 * self.e_mono
                self.stats["evaluated"] += int(live.sum())
            if self.charges is not None and self.cutoff is not None:
                qq = self.charges[:, None] * self.charges[None, :]
                w = _switch(d, self.cutoff - self.lr_window, self.cutoff)
                lr = (KE * qq * w / np.maximum(d, 1e-9)).sum(axis=(1, 2))
                chunk_out[live] += lr[live]
            self.stats["walled"] += int(walled.sum())
            out[start:stop] = chunk_out
            if self.progress is not None:
                self.progress(r, stop, total)
        return out

    def slab(self, r, ang):
        na = len(ang)
        full = na - 1
        if full % self.order:
            raise ValueError(f"C{self.order} does not divide the {full}-step angle grid")
        p = full // self.order
        uniq, inverse = self._orbit(p)
        angles = np.deg2rad(np.arange(p) * self.da)
        t1 = angles[uniq // p]
        t2 = angles[uniq % p]
        values = self._energies(r, t1, t2)
        torus = values[inverse]
        idx = np.arange(na) % p
        self.stats["rows"] += 1
        return torus[np.ix_(idx, idx)]
