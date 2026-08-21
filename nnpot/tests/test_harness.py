import tempfile
import unittest
from pathlib import Path

import numpy as np

import nnpot
from fsmp_gui.generate import GridSpec, MMFFBackend, generate
from fsmp_gui.molecule import Atom, Molecule
from nnpot.backend import CAP_JMOL, KE, NNBackend, _switch
from nnpot.geometry import (ATOMIC_NUMBER, dimer_coords, orbit_map,
                           rotational_order, symmetrize, symmetry_error)
from nnpot.reference import MMFFCalculator

MODEL = Path(__file__).resolve().parents[2] / "samples" / "models" / "trimesic_acid.xyz"
DA = 15.0


def angles(da=DA):
    return np.deg2rad(np.arange(int(round(360.0 / da)) + 1) * da)


def tma():
    return Molecule.load_xyz(MODEL)


def c3_tma():
    mol = tma()
    elements = [a.element for a in mol.atoms]
    xy = symmetrize(np.array([(a.x, a.y) for a in mol.atoms]), elements, 3)
    return Molecule([Atom(e, float(xy[i, 0]), float(xy[i, 1]), 0.0)
                     for i, e in enumerate(elements)], mol.comment)


def harness(mol, order, wall=0.0):
    return NNBackend([a.element for a in mol.atoms],
                     np.array([(a.x, a.y) for a in mol.atoms]),
                     MMFFCalculator(mol), DA, order=order, wall=wall)


class HarnessTest(unittest.TestCase):

    def test_reproduces_the_mmff_backend(self):
        mol = tma()
        reference = MMFFBackend(mol)
        backend = harness(mol, order=1)
        for r in (8.0, 10.0, 13.0):
            with self.subTest(r=r):
                self.assertTrue(np.allclose(backend.slab(r, angles()),
                                            reference.slab(r, angles()),
                                            rtol=1e-9, atol=1e-6))

    def test_symmetrized_molecule_is_exactly_c3(self):
        mol = tma()
        elements = [a.element for a in mol.atoms]
        raw = np.array([(a.x, a.y) for a in mol.atoms])
        self.assertGreater(symmetry_error(raw, elements, 3), 0.0)
        sym = np.array([(a.x, a.y) for a in c3_tma().atoms])
        self.assertLess(symmetry_error(sym, elements, 3), 1e-9)
        self.assertEqual(rotational_order(sym, elements), 3)

    def test_c3_folding_gives_the_same_grid(self):
        mol = c3_tma()
        reference = MMFFBackend(mol)
        backend = harness(mol, order=3)
        for r in (8.0, 11.0):
            with self.subTest(r=r):
                self.assertTrue(np.allclose(backend.slab(r, angles()),
                                            reference.slab(r, angles()),
                                            rtol=1e-9, atol=1e-6))

    def test_exchange_symmetry_is_exact(self):
        slab = harness(c3_tma(), order=3).slab(9.0, angles())
        na = slab.shape[0]
        full = na - 1
        idx = (np.arange(na) + full // 2) % full
        partner = slab[np.ix_(idx, idx)].T
        self.assertTrue(np.allclose(slab, partner, rtol=1e-9, atol=1e-6))

    def test_wall_caps_overlapping_configurations(self):
        backend = harness(c3_tma(), order=3, wall=3.0)
        slab = backend.slab(6.0, angles())
        self.assertTrue(np.all(slab == CAP_JMOL))
        self.assertEqual(backend.stats["evaluated"], 0)

    def test_orbit_map_halves_the_configurations(self):
        uniq, inverse = orbit_map(24, DA)
        self.assertEqual(inverse.shape, (24, 24))
        self.assertLess(len(uniq), 0.55 * 24 * 24)
        self.assertEqual(set(np.unique(inverse)), set(range(len(uniq))))

    def test_generated_file_matches_the_mmff_generator(self):
        mol = c3_tma()
        spec = GridSpec(r_min=8.0, r_max=9.0, dr=0.5, da=DA, fold_deg=120.0)
        with tempfile.TemporaryDirectory() as tmp:
            mine, theirs = Path(tmp) / "nn.bin", Path(tmp) / "mmff.bin"
            self.assertTrue(generate(harness(mol, order=3), spec, mine))
            self.assertTrue(generate(MMFFBackend(mol), spec, theirs))
            self.assertEqual(mine.read_bytes(), theirs.read_bytes())


class _Coulomb:

    def __init__(self, q, cutoff=None, window=1.0):
        self.q = np.asarray(q, dtype=float)
        self.n = len(self.q)
        self.cutoff = cutoff
        self.window = window

    def energies(self, z, coords):
        coords = np.asarray(coords, dtype=float)
        k, m, _ = coords.shape
        if m == self.n:
            return np.zeros(k)
        a, b = coords[:, :self.n, :], coords[:, self.n:, :]
        delta = a[:, :, None, :] - b[:, None, :, :]
        d = np.maximum(np.sqrt(np.einsum("kijc,kijc->kij", delta, delta)), 1e-9)
        qq = (self.q[:, None] * self.q[None, :])[None, :, :]
        weight = 1.0
        if self.cutoff is not None:
            weight = 1.0 - _switch(d, self.cutoff - self.window, self.cutoff)
        return (KE * qq * weight / d).sum(axis=(1, 2))


class LongRangeCorrectionTest(unittest.TestCase):

    def setUp(self):
        from fsmp_gui.mmff import mmff_pair_params
        self.mol = c3_tma()
        self.elements = [a.element for a in self.mol.atoms]
        self.xy = np.array([(a.x, a.y) for a in self.mol.atoms])
        self.q = np.asarray(mmff_pair_params(self.mol).q, dtype=float)

    def backend(self, calculator, charges):
        return NNBackend(self.elements, self.xy, calculator, DA,
                         charges=charges, order=3, wall=0.0, lr_window=1.0)

    def test_truncated_model_plus_correction_is_the_complete_model(self):
        cut = self.backend(_Coulomb(self.q, cutoff=5.0), self.q)
        whole = self.backend(_Coulomb(self.q), None)
        for r in (8.0, 12.0, 20.0):
            with self.subTest(r=r):
                self.assertTrue(np.allclose(cut.slab(r, angles()),
                                            whole.slab(r, angles()),
                                            rtol=1e-9, atol=1e-6))
        self.assertGreater(cut.stats["rows_beyond_reach"], 0)

    def test_a_model_without_a_cutoff_is_left_alone(self):
        whole = self.backend(_Coulomb(self.q), self.q)
        bare = self.backend(_Coulomb(self.q), None)
        self.assertTrue(np.allclose(whole.slab(9.0, angles()),
                                    bare.slab(9.0, angles()),
                                    rtol=1e-9, atol=1e-6))
        self.assertEqual(whole.stats["rows_beyond_reach"], 0)


class AIMNet2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        try:
            import aimnet2calc  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("aimnet2calc is not installed")
        from nnpot.aimnet import AIMNet2
        cls.calc = AIMNet2()
        cls.mol = c3_tma()
        cls.xy = np.array([(a.x, a.y) for a in cls.mol.atoms])
        cls.z = np.array([ATOMIC_NUMBER[a.element] for a in cls.mol.atoms])

    def test_simple_coulomb_means_no_declared_cutoff(self):
        described = self.calc.describe()
        self.assertEqual(described["coulomb"], "simple")
        self.assertIsNone(described["declared_cutoff"])
        self.assertGreater(described["short_cutoff"], 0.0)

    def test_far_apart_molecules_barely_interact(self):
        mono = np.concatenate([self.xy, np.zeros((len(self.xy), 1))], axis=1)
        alone = self.calc.energies(self.z, mono[None])[0]
        t = np.deg2rad(np.array([0.0, 37.0, 91.0, 150.0]))
        coords = dimer_coords(self.xy, 60.0, t, t[::-1])
        pair = self.calc.energies(np.concatenate([self.z, self.z]), coords)
        self.assertLess(np.abs(pair - 2.0 * alone).max(), 20.0)


if __name__ == "__main__":
    unittest.main()
