"""Tests for the non-GUI core: xyz input/output and the project manifest.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_core.py
"""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.forcefield import ForcefieldError, read_header
from fsmp_gui.molecule import Atom, Molecule
from fsmp_gui.project import Project, ProjectError, safe_filename
from fsmp_gui.sitemodel import Site, SiteModel

REPO = Path(__file__).resolve().parents[2]
TEST_GRID = REPO / "tests" / "data" / "TMA_simple_2020_s4.v2.bin"


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


class TestSiteModel(unittest.TestCase):
    def test_roundtrip(self):
        model = SiteModel([Site("O", 0.0, 0.0, 0.0, -0.5, 500.0, 3.1),
                           Site("q+", 1.2, 0.0, 0.0, 0.5, 0.0, 0.0)],
                          "test model")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "m.site"
            model.save(path)
            back = SiteModel.load(path)
        self.assertEqual(len(back.sites), 2)
        self.assertEqual(back.sites[0].label, "O")
        self.assertAlmostEqual(back.sites[0].epsilon, 500.0)
        self.assertAlmostEqual(back.total_charge(), 0.0)
        self.assertTrue(back.sites[0].is_lj)
        self.assertFalse(back.sites[1].is_lj)
        self.assertTrue(back.sites[1].is_charge)

    def test_load_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.site"
            bad.write_text("1\ncomment\nO 0 0 0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                SiteModel.load(bad)


class TestProject(unittest.TestCase):
    def test_create_open_set_atomistic(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "proj"
            project = Project.create(root, "test project")
            self.assertIsNone(project.model)
            project.set_atomistic("benzene", Molecule([Atom("C", 0, 0)]))

            again = Project.open(root)
            self.assertEqual(again.name, "test project")
            self.assertEqual(again.model["name"], "benzene")
            self.assertEqual(again.model_kind, "atomistic")
            self.assertTrue(again.model_path().is_file())

            again.clear_model()
            self.assertIsNone(again.model)
            self.assertFalse((root / "molecules" / "benzene.xyz").exists())

    def test_set_site_model(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_site_model("sites", SiteModel([Site("LJ", 0, 0)]))
            self.assertEqual(project.model_kind, "site")
            self.assertTrue(project.model_path().name.endswith(".site"))
            again = Project.open(project.root)
            self.assertEqual(again.model_kind, "site")

    def test_switch_model_kind_deletes_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_atomistic("first", Molecule([Atom("C", 0, 0)]))
            old_path = project.model_path()
            project.set_site_model("second", SiteModel([Site("LJ", 0, 0)]))
            self.assertFalse(old_path.exists())
            self.assertEqual(project.model_kind, "site")
            self.assertTrue(project.model_path().is_file())

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

    def test_charges_lifecycle(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_atomistic("m", Molecule([Atom("C", 0, 0), Atom("O", 1, 0)]))
            project.set_charges("Gasteiger", [0.1, -0.1])
            again = Project.open(project.root)
            self.assertEqual(again.charges["method"], "Gasteiger")
            self.assertEqual(again.charges["values"], [0.1, -0.1])
            # replacing the molecule invalidates the stored charges
            again.set_atomistic("m", Molecule([Atom("C", 0, 0)]))
            self.assertIsNone(again.charges)

    def test_potential_slot(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            self.assertIsNone(project.potential)

            outside = Path(tmp) / "ff.bin"
            outside.write_bytes(b"x")
            project.set_potential("ff", outside)
            self.assertTrue(Path(project.potential["path"]).is_absolute())
            self.assertEqual(project.potential_path(), outside.resolve())

            inside = project.root / "local.bin"
            inside.write_bytes(b"x")
            project.set_potential("local", inside)
            self.assertEqual(project.potential["path"], "local.bin")
            self.assertTrue(project.potential_path().is_file())

            project.clear_potential()
            self.assertIsNone(project.potential)
            self.assertTrue(inside.is_file())  # detaching never deletes


class TestForcefieldHeader(unittest.TestCase):
    @unittest.skipUnless(TEST_GRID.is_file(), "committed test grid not found")
    def test_reads_committed_grid(self):
        info = read_header(TEST_GRID)
        self.assertEqual(info.version, 2)
        self.assertEqual(info.dtype, "float")
        self.assertEqual(info.n_dist, 281)
        self.assertEqual(info.n_ang, 31)
        self.assertAlmostEqual(info.min_dist, 7.6)
        self.assertAlmostEqual(info.fold, 120.0)
        self.assertIn("float", info.summary())

    def test_rejects_garbage(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.bin"
            bad.write_bytes(b"NOPE" + b"\0" * 60)
            with self.assertRaises(ForcefieldError):
                read_header(bad)
            short = Path(tmp) / "short.bin"
            short.write_bytes(b"FSMP")
            with self.assertRaises(ForcefieldError):
                read_header(short)


if __name__ == "__main__":
    unittest.main(verbosity=2)
