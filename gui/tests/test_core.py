"""Tests for the non-GUI core: xyz input/output and the project manifest.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_core.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.molecule import Atom, Molecule
from fsmp_gui.project import Project, ProjectError, safe_filename


class TestMolecule(unittest.TestCase):
    def test_xyz_roundtrip(self):
        mol = Molecule([Atom("C", 0.0, 0.0, 0.0), Atom("O", 1.16, 0.0, 0.0),
                        Atom("O", -1.16, 0.0, 0.0)], "carbon dioxide")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "co2.xyz"
            mol.save_xyz(path)
            back = Molecule.load_xyz(path)
        self.assertEqual(len(back.atoms), 3)
        self.assertEqual(back.atoms[1].element, "O")
        self.assertAlmostEqual(back.atoms[1].x, 1.16)
        self.assertEqual(back.comment, "carbon dioxide")

    def test_extra_columns_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "m.xyz"
            path.write_text("1\nwith charge column\nN 1.0 2.0 3.0 -0.8\n",
                            encoding="utf-8")
            mol = Molecule.load_xyz(path)
        self.assertEqual(mol.atoms[0].element, "N")
        self.assertAlmostEqual(mol.atoms[0].z, 3.0)

    def test_formula_hill_order(self):
        mol = Molecule([Atom("O", 0, 0), Atom("H", 0, 0), Atom("C", 0, 0),
                        Atom("H", 0, 0), Atom("O", 0, 0)])
        self.assertEqual(mol.formula(), "CH2O2")

    def test_load_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.xyz"
            bad.write_text("2\ncomment\nC 0 0 0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                Molecule.load_xyz(bad)


class TestProject(unittest.TestCase):
    def test_create_open_set_molecule(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "proj"
            project = Project.create(root, "test project")
            self.assertIsNone(project.molecule)
            project.set_molecule("benzene", Molecule([Atom("C", 0, 0)]))

            again = Project.open(root)
            self.assertEqual(again.name, "test project")
            self.assertEqual(again.molecule["name"], "benzene")
            self.assertTrue(again.molecule_path().is_file())

            again.clear_molecule()
            self.assertIsNone(again.molecule)
            self.assertFalse((root / "molecules" / "benzene.xyz").exists())

    def test_replace_molecule_deletes_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_molecule("first", Molecule([Atom("C", 0, 0)]))
            old_path = project.molecule_path()
            project.set_molecule("second", Molecule([Atom("O", 0, 0)]))
            self.assertFalse(old_path.exists())
            self.assertEqual(project.molecule["name"], "second")
            self.assertTrue(project.molecule_path().is_file())

    def test_open_rejects_non_project(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ProjectError):
                Project.open(tmp)

    def test_create_rejects_non_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "something.txt").write_text("x", encoding="utf-8")
            with self.assertRaises(ProjectError):
                Project.create(tmp, "x")

    def test_safe_filename(self):
        self.assertEqual(safe_filename('a b/c:d*e'), "a_b_c_d_e")
        self.assertEqual(safe_filename("  "), "unnamed")


if __name__ == "__main__":
    unittest.main(verbosity=2)
