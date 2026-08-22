import tempfile
import unittest
from pathlib import Path

import numpy as np

import nnpot
from fsmp_gui.generate import GridSpec, MMFFBackend, generate
from fsmp_gui.molecule import Atom, Molecule, rotated
from nnpot.compare import ComparisonError, compare
from nnpot.geometry import symmetrize

MODEL = Path(__file__).resolve().parents[2] / "samples" / "models" / "trimesic_acid.xyz"
DA = 15.0
SPEC = GridSpec(r_min=8.0, r_max=8.4, dr=0.2, da=DA, fold_deg=120.0)


def c3_tma():
    mol = Molecule.load_xyz(MODEL)
    elements = [a.element for a in mol.atoms]
    xy = symmetrize(np.array([(a.x, a.y) for a in mol.atoms]), elements, 3)
    return Molecule([Atom(e, float(xy[i, 0]), float(xy[i, 1]), 0.0)
                     for i, e in enumerate(elements)], mol.comment)


class CompareTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.molecule = c3_tma()

    def write(self, name, molecule, spec=SPEC):
        path = Path(self.tmp.name) / name
        self.assertTrue(generate(MMFFBackend(molecule), spec, path))
        return path

    def test_a_grid_matches_itself(self):
        path = self.write("one.bin", self.molecule)
        result = compare(path, path)
        self.assertEqual(result["shift"], 0.0)
        self.assertFalse(result["mirror"])
        self.assertLess(result["rmse"], 1e-6)

    def test_a_turned_molecule_is_recognised(self):
        straight = self.write("straight.bin", self.molecule)
        turned = self.write("turned.bin",
                            Molecule(rotated(self.molecule.atoms, 30.0)))
        result = compare(straight, turned)
        self.assertLess(result["rmse"], 1.0)
        self.assertAlmostEqual(result["shift"] % 120.0, 90.0, places=6)

    def test_mismatched_folding_is_refused(self):
        folded = self.write("folded.bin", self.molecule)
        whole = self.write("whole.bin", self.molecule,
                           GridSpec(r_min=8.0, r_max=8.4, dr=0.2, da=DA,
                                    fold_deg=360.0))
        with self.assertRaises(ComparisonError):
            compare(folded, whole)


if __name__ == "__main__":
    unittest.main()
