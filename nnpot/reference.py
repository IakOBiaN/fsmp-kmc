import numpy as np

KCAL_TO_JMOL = 4184.0
MMFF_ELE = 332.0716
MMFF_DELTA = 0.05


class MMFFCalculator:

    name = "MMFF94"
    cutoff = None

    def __init__(self, molecule, net_charge=0):
        from fsmp_gui.mmff import mmff_pair_params
        params = mmff_pair_params(molecule, net_charge)
        self.q = np.asarray(params.q, dtype=float)
        self.Rstar = np.asarray(params.Rstar, dtype=float)
        self.eps = np.asarray(params.eps, dtype=float)
        self.n = len(self.q)

    def energies(self, z, coords):
        coords = np.asarray(coords, dtype=float)
        k, m, _ = coords.shape
        if m == self.n:
            return np.zeros(k)
        if m != 2 * self.n:
            raise ValueError("MMFFCalculator scores one or two copies of its molecule")
        a = coords[:, :self.n, :]
        b = coords[:, self.n:, :]
        delta = a[:, :, None, :] - b[:, None, :, :]
        d = np.sqrt(np.einsum("kijc,kijc->kij", delta, delta))
        np.maximum(d, 1e-6, out=d)
        rs = self.Rstar[None, :, :]
        t = (1.07 * rs) / (d + 0.07 * rs)
        t7 = t ** 7
        r7 = rs ** 7
        vdw = self.eps[None, :, :] * t7 * (1.12 * r7 / (d ** 7 + 0.12 * r7) - 2.0)
        qq = (self.q[:, None] * self.q[None, :])[None, :, :]
        ele = MMFF_ELE * qq / (d + MMFF_DELTA)
        return (vdw + ele).sum(axis=(1, 2)) * KCAL_TO_JMOL
