"""Tests for the RDKit-backed chemistry layer (skipped without RDKit).

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_chem.py
"""

import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui import chem
from fsmp_gui.molecule import Atom, Molecule


def benzene() -> Molecule:
    atoms = []
    for i in range(6):
        a = math.radians(60 * i)
        atoms.append(Atom("C", 1.42 * math.cos(a), 1.42 * math.sin(a)))
        atoms.append(Atom("H", 2.60 * math.cos(a), 2.60 * math.sin(a)))
    return Molecule(atoms, "benzene, slightly stretched")


@unittest.skipUnless(chem.HAVE_RDKIT, "RDKit not installed")
class TestChem(unittest.TestCase):
    def test_optimize_lowers_energy(self):
        res = chem.optimize(benzene(), "MMFF94")
        self.assertTrue(res.converged)
        self.assertLess(res.e_after, res.e_before)
        self.assertEqual(len(res.molecule.atoms), 12)
        self.assertEqual(res.molecule.atoms[0].element, "C")
        # C-C bond should relax towards the aromatic 1.39 A
        a, b = res.molecule.atoms[0], res.molecule.atoms[2]
        d = math.dist((a.x, a.y, a.z), (b.x, b.y, b.z))
        self.assertAlmostEqual(d, 1.39, delta=0.03)

    def test_optimize_uff(self):
        res = chem.optimize(benzene(), "UFF")
        self.assertLessEqual(res.e_after, res.e_before)

    def test_gasteiger_charges(self):
        values = chem.partial_charges(benzene(), "Gasteiger")
        self.assertEqual(len(values), 12)
        self.assertAlmostEqual(sum(values), 0.0, places=3)
        # hydrogens positive, carbons negative in benzene
        self.assertGreater(values[1], 0)
        self.assertLess(values[0], 0)

    def test_mmff_charges(self):
        values = chem.partial_charges(benzene(), "MMFF94")
        self.assertEqual(len(values), 12)
        self.assertAlmostEqual(sum(values), 0.0, places=3)

    def test_empty_molecule_error(self):
        with self.assertRaises(chem.ChemError):
            chem.optimize(Molecule())

    def test_unknown_element_error(self):
        with self.assertRaises(chem.ChemError):
            chem.optimize(Molecule([Atom("Qq", 0, 0)]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
