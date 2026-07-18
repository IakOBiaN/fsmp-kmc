"""The RDKit backend: optimizing a geometry and generating an MMFF94 pair
potential. The key physical guarantee is that MMFF94 finds a deep, attractive
hydrogen-bond well at the right geometry, so a bonded assembly does not fly
apart. Skipped entirely when RDKit is not installed.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_mmff.py
"""

import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from fsmp_gui import selftest
from fsmp_gui.engine import find_engine
from fsmp_gui.forcefield import ForcefieldGrid, read_header
from fsmp_gui.generate import GridSpec, MMFFBackend, generate
from fsmp_gui.mmff import (MMFFError, _rdkit_mol, mmff_pair_params,
                           optimize_molecule, rdkit_available)
from fsmp_gui.molecule import Atom, Molecule, bonded_pairs, connected_components

REPO = Path(__file__).resolve().parents[2]
TMA = REPO / "models" / "trimesic_acid.xyz"

skip_no_rdkit = unittest.skipUnless(rdkit_available(), "RDKit not installed")

# a compressed trimesic-acid geometry (short ring bonds): the ring is squished
# enough that RDKit's own perceiver invents spurious cross-ring bonds and fails
_DISTORTED_TMA = [
    ("C", -1.40817834, 2.51673644), ("O", -0.87542713, 3.61018695),
    ("O", -2.74488711, 2.38625041), ("C", -0.70511375, 1.21103238),
    ("C", -1.02570000, 0.05420000), ("C", -0.69622801, -1.21616310),
    ("C", 0.70410080, -1.20980425), ("C", 0.85010000, 0.04450000),
    ("C", 0.69567047, 1.21467066), ("H", -1.97710000, -0.48880000),
    ("H", 1.23037969, 2.16365787), ("H", 1.25859345, -2.14736948),
    ("C", -1.47546849, -2.47788697), ("O", -2.68880001, -2.56323597),
    ("O", -0.69410985, -3.57026773), ("O", 3.56422807, -1.04695054),
    ("C", 2.50960000, -0.01920000), ("O", 3.14370000, 1.34150000),
    ("H", -3.07823289, 3.30848281), ("H", 3.61690000, 2.23220000),
    ("H", -1.32611373, -4.32006980)]


def distorted_tma() -> Molecule:
    return Molecule([Atom(el, x, y, 0.0) for el, x, y in _DISTORTED_TMA])


def _coarse_well(backend, da=5.0, r_lo=6.0, r_hi=14.0, dr=0.5):
    """Global minimum (kJ/mol) and its (r, a1, a2) over a coarse scan. The
    grid is deliberately coarse: it finds the same TMA well (-42 kJ/mol at
    r = 10) as a fine scan, in a fraction of the time."""
    ang = np.deg2rad(np.arange(int(round(360 / da)) + 1) * da)
    best = (1e18, None)
    for r in np.arange(r_lo, r_hi + 1e-9, dr):
        slab = backend.slab(r, ang)
        k = np.unravel_index(int(np.argmin(slab)), slab.shape)
        if slab[k] < best[0]:
            best = (float(slab[k]), (r, math.degrees(ang[k[0]]),
                                     math.degrees(ang[k[1]])))
    return best[0] / 1000.0, best[1]


@skip_no_rdkit
class TestOptimize(unittest.TestCase):
    def test_tma_relaxes_planar_and_centred(self):
        mol = Molecule.load_xyz(TMA)
        relaxed, report = optimize_molecule(mol)

        self.assertEqual(len(relaxed.atoms), len(mol.atoms))
        self.assertEqual([a.element for a in relaxed.atoms],
                         [a.element for a in mol.atoms])
        self.assertEqual(report.method, "MMFF94")
        self.assertTrue(report.converged)

        # kept planar and recentred on the molecule centre (the engine frame)
        self.assertTrue(all(a.z == 0.0 for a in relaxed.atoms))
        self.assertLess(report.max_out_of_plane, 0.1)   # TMA is flat
        cx = sum(a.x for a in relaxed.atoms) / len(relaxed.atoms)
        cy = sum(a.y for a in relaxed.atoms) / len(relaxed.atoms)
        self.assertAlmostEqual(cx, 0.0, places=6)
        self.assertAlmostEqual(cy, 0.0, places=6)

        # geometry stays sane: no collapsed atoms, comparable overall size
        pts = np.array([(a.x, a.y) for a in relaxed.atoms])
        d = np.hypot(pts[:, None, 0] - pts[None, :, 0],
                     pts[:, None, 1] - pts[None, :, 1])
        d[np.diag_indices_from(d)] = np.inf
        self.assertGreater(d.min(), 0.9)                # closest atoms > 0.9 A
        r_ext = np.hypot(pts[:, 0], pts[:, 1]).max()
        self.assertLess(r_ext, 6.0)                     # TMA radius ~4.7 A

    def test_distorted_geometry_is_perceived(self):
        """A compressed trimesic-acid geometry (short ring bonds) is what the
        user hit: it must still perceive and optimize into a proper molecule."""
        relaxed, report = optimize_molecule(distorted_tma())
        self.assertEqual(len(relaxed.atoms), 21)
        self.assertEqual(report.method, "MMFF94")
        # the relaxed molecule is a proper, MMFF-typable trimesic acid
        self.assertEqual(mmff_pair_params(relaxed).xy.shape, (21, 2))

    def test_editor_bonds_match_perception(self):
        """The bonds the canvas draws (bonded_pairs) are exactly the connectivity
        the optimizer perceives, for a clean and a compressed geometry. This is
        the what-you-see-is-what-you-get guarantee."""
        for mol in (Molecule.load_xyz(TMA), distorted_tma()):
            drawn = set(bonded_pairs(mol.atoms))
            perceived = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx())))
                         for b in _rdkit_mol(mol, 0).GetBonds()}
            self.assertEqual(drawn, perceived, msg=f"{len(mol.atoms)} atoms")

    def test_disconnected_molecule_is_rejected(self):
        """A molecule model must be whole: an atom detached from the rest is
        flagged, not silently bridged."""
        atoms = [Atom(el, x, y, 0.0) for el, x, y in _DISTORTED_TMA]
        atoms[-1] = Atom("H", 25.0, 25.0, 0.0)   # a stray hydrogen, far away
        mol = Molecule(atoms)
        self.assertEqual(len(connected_components(mol.atoms)), 2)
        with self.assertRaises(MMFFError):
            optimize_molecule(mol)

    def test_empty_molecule_is_rejected(self):
        with self.assertRaises(MMFFError):
            optimize_molecule(Molecule())

    def test_garbage_geometry_is_rejected(self):
        # three atoms piled on top of each other: bond perception must fail
        # cleanly rather than crash
        blob = Molecule([Atom("C", 0.0, 0.0, 0.0), Atom("C", 0.01, 0.0, 0.0),
                         Atom("C", 0.0, 0.01, 0.0)])
        with self.assertRaises(MMFFError):
            optimize_molecule(blob)


@skip_no_rdkit
class TestPairEnergy(unittest.TestCase):
    def test_typed_parameters(self):
        params = mmff_pair_params(Molecule.load_xyz(TMA))
        n = len(params.elements)
        self.assertEqual(params.xy.shape, (n, 2))
        self.assertEqual(params.Rstar.shape, (n, n))
        # trimesic acid is neutral: MMFF partial charges sum to ~0
        self.assertAlmostEqual(float(params.q.sum()), 0.0, places=4)

    def test_hydrogen_bond_well_is_deep(self):
        """The whole point: MMFF94 must produce a deep attractive well near the
        carboxyl dimer geometry (validated ~-43 kJ/mol at r~10). This is what
        keeps a hydrogen-bonded crystal from disintegrating."""
        backend = MMFFBackend(Molecule.load_xyz(TMA))
        well, (r, a1, a2) = _coarse_well(backend)
        self.assertLess(well, -30.0)          # firmly bound (>> kT)
        self.assertGreater(well, -60.0)        # but not unphysically deep
        self.assertTrue(9.0 <= r <= 11.0, msg=f"well at r={r}")

    def test_tails_and_core(self):
        backend = MMFFBackend(Molecule.load_xyz(TMA))
        ang = np.deg2rad(np.arange(0, 361, 30))
        far = backend.slab(30.0, ang)          # dispersion tail, small but real
        self.assertLess(np.abs(far).max() / 1000.0, 1.0)
        near = backend.slab(4.0, ang)          # deep overlap: strongly repulsive
        self.assertGreater(near.max() / 1000.0, 100.0)

    def test_net_charge_species(self):
        """The generator's net-charge path: a planar acetate anion perceives
        with charge -1, MMFF types it, and its partial charges sum to -1."""
        acetate = Molecule([Atom(el, x, y, 0.0) for el, x, y in
                            [("C", -1.50, 0.00), ("C", 0.00, 0.00),
                             ("O", 0.65, 1.05), ("O", 0.65, -1.05),
                             ("H", -2.10, 0.90), ("H", -2.10, -0.90),
                             ("H", -2.59, 0.00)]])
        params = mmff_pair_params(acetate, net_charge=-1)
        self.assertAlmostEqual(float(params.q.sum()), -1.0, places=4)
        backend = MMFFBackend(acetate, net_charge=-1)
        slab = backend.slab(8.0, np.deg2rad(np.arange(0, 361, 30)))
        self.assertTrue(np.isfinite(slab).all())
        self.assertGreater(float(slab.min()), 0.0)   # two anions repel


@skip_no_rdkit
class TestSelftest(unittest.TestCase):
    @unittest.skipUnless(find_engine(), "engine binary not present")
    def test_reports_ok_from_source(self):
        """The bundle health check must pass from a source tree too: it is
        the same code path the release workflow runs on every frozen bundle."""
        with tempfile.TemporaryDirectory() as td:
            report = Path(td) / "selftest.log"
            with contextlib.redirect_stdout(io.StringIO()):
                code = selftest.run(str(report))
            self.assertEqual(code, 0)
            text = report.read_text(encoding="utf-8")
        self.assertIn("SELFTEST OK", text)
        self.assertIn("rdkit", text)
        self.assertIn("engine", text)


@skip_no_rdkit
class TestGenerateAtomistic(unittest.TestCase):
    def test_round_trip_v2_file(self):
        backend = MMFFBackend(Molecule.load_xyz(TMA))
        spec = GridSpec(6.0, 13.0, 0.2, 5.0, 120.0, use_float=True)
        # ignore_cleanup_errors: ForcefieldGrid memory-maps the file, so on
        # Windows it may still be mapped when the temp dir is removed
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
            out = Path(td) / "tma_mmff.v2.bin"
            generate(backend, spec, out)
            info = read_header(out)
            self.assertEqual(info.version, 2)
            self.assertEqual(info.n_dist, spec.n_dist)
            self.assertEqual(info.n_ang, int(round(120 / spec.da)) + 1)
            self.assertAlmostEqual(info.fold, 120.0)

            grid = ForcefieldGrid.open(out)
            self.assertLess(float(np.asarray(grid._data).min()) / 1000.0, -25.0)
            # beyond the cutoff the reader returns exactly zero
            self.assertEqual(grid.energy_at(info.r_max + 5.0, 0, 0), 0.0)
            grid._data._mmap.close()   # release the map before cleanup


if __name__ == "__main__":
    unittest.main(verbosity=2)
