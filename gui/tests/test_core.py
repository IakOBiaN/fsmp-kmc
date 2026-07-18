"""Tests for the non-GUI core: xyz input/output and the project manifest.

Run from the repository root:
    gui/.venv/Scripts/python gui/tests/test_core.py
"""

import math
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsmp_gui.forcefield import ForcefieldError, read_header
from fsmp_gui.molecule import (Atom, Molecule, aimed_at_x, centroid, rotated,
                               translated)
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


class TestGeometry(unittest.TestCase):
    """Whole-molecule transforms behind the editor's Center / rotate /
    Point-at-+x actions: rigid, z-preserving, engine-frame conventions."""

    ATOMS = [Atom("C", 1.2, 0.0, 0.3), Atom("O", 0.0, 3.4, -1.0),
             Atom("H", -2.0, 0.7)]

    def assert_rigid(self, a, b):
        """Same elements, same pairwise (3D) distances."""
        self.assertEqual([x.element for x in a], [x.element for x in b])
        for i in range(len(a)):
            for j in range(i + 1, len(a)):
                da = math.dist((a[i].x, a[i].y, a[i].z),
                               (a[j].x, a[j].y, a[j].z))
                db = math.dist((b[i].x, b[i].y, b[i].z),
                               (b[j].x, b[j].y, b[j].z))
                self.assertAlmostEqual(da, db, places=12, msg=f"pair {i},{j}")

    def test_centroid_and_translate(self):
        atoms = [Atom("C", 1.0, 2.0, 0.5), Atom("O", 3.0, 6.0)]
        cx, cy = centroid(atoms)
        self.assertEqual((cx, cy), (2.0, 4.0))
        moved = translated(atoms, -cx, -cy)
        self.assertEqual(centroid(moved), (0.0, 0.0))
        self.assertEqual(moved[0].z, 0.5)                  # z rides along
        self.assert_rigid(atoms, moved)

    def test_rotation_is_rigid_and_counterclockwise(self):
        turned = rotated(self.ATOMS, 37.0)
        self.assert_rigid(self.ATOMS, turned)
        self.assertEqual(turned[1].z, -1.0)
        # a positive angle turns +x toward +y
        e = rotated([Atom("C", 1.0, 0.0)], 90.0)[0]
        self.assertAlmostEqual(e.x, 0.0, places=12)
        self.assertAlmostEqual(e.y, 1.0, places=12)

    def test_aim_at_x_lands_exactly_on_the_axis(self):
        for index in range(len(self.ATOMS)):
            aimed = aimed_at_x(self.ATOMS, index)
            a = self.ATOMS[index]
            self.assertEqual(aimed[index].y, 0.0)
            self.assertAlmostEqual(aimed[index].x, math.hypot(a.x, a.y),
                                   places=12)
            self.assertGreater(aimed[index].x, 0.0)
            self.assert_rigid(self.ATOMS, aimed)

    def test_aim_at_x_refuses_the_origin(self):
        with self.assertRaises(ValueError):
            aimed_at_x([Atom("C", 0.0, 0.0), Atom("O", 1.0, 0.0)], 0)


class TestSiteModel(unittest.TestCase):
    def test_roundtrip(self):
        model = SiteModel([Site("O", 0.0, 0.0, 0.0, -0.5, 500.0, 3.1, 2.4),
                           Site("q+", 1.2, 0.0, 0.0, 0.5, 0.0, 0.0)],
                          "test model")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "m.site"
            model.save(path)
            back = SiteModel.load(path)
        self.assertEqual(len(back.sites), 2)
        self.assertEqual(back.sites[0].label, "O")
        self.assertAlmostEqual(back.sites[0].epsilon, 500.0)
        self.assertAlmostEqual(back.sites[0].r0, 2.4)
        self.assertAlmostEqual(back.total_charge(), 0.0)
        self.assertTrue(back.sites[0].is_lj)
        self.assertFalse(back.sites[1].is_lj)
        self.assertTrue(back.sites[1].is_charge)

    def test_load_without_r0_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "old.site"
            path.write_text("1\nold 7-column file\nO 0 0 0 -0.5 500 3.1\n",
                            encoding="utf-8")
            back = SiteModel.load(path)
        self.assertAlmostEqual(back.sites[0].sigma, 3.1)
        self.assertAlmostEqual(back.sites[0].r0, 0.0)

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
            self.assertIsNone(project.atomistic)
            project.set_atomistic("benzene", Molecule([Atom("C", 0, 0)]))

            again = Project.open(root)
            self.assertEqual(again.name, "test project")
            self.assertEqual(again.atomistic["name"], "benzene")
            self.assertEqual(again.generation_kind, "atomistic")
            self.assertTrue(again.model_path(again.atomistic).is_file())

            again.clear_atomistic()
            self.assertIsNone(again.atomistic)
            self.assertFalse((root / "molecules" / "benzene.xyz").exists())

    def test_two_models_coexist(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_atomistic("benzene", Molecule([Atom("C", 0, 0)]))
            project.set_site("sites", SiteModel([Site("LJ", 0, 0)]))
            again = Project.open(project.root)
            # both slots are independent and both are kept
            self.assertEqual(again.atomistic["name"], "benzene")
            self.assertEqual(again.site["name"], "sites")
            self.assertEqual(again.generation_kind, "site")   # site wins
            self.assertTrue(again.model_path(again.atomistic).is_file())
            self.assertTrue(again.model_path(again.site).name.endswith(".site"))
            # clearing one leaves the other
            again.clear_site()
            self.assertIsNone(again.site)
            self.assertIsNotNone(again.atomistic)
            self.assertEqual(again.generation_kind, "atomistic")

    def test_replace_within_slot_deletes_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            project.set_site("first", SiteModel([Site("LJ", 0, 0)]))
            old_path = project.model_path(project.site)
            project.set_site("second", SiteModel([Site("LJ", 0, 0)]))
            self.assertFalse(old_path.exists())
            self.assertEqual(project.site["name"], "second")

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

    def test_unit_cell_slot(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Project.create(Path(tmp) / "proj", "p")
            self.assertIsNone(project.unit_cell)
            project.set_unit_cell(11.1, 19.2, [{"x": 0, "y": 0, "phi": 0},
                                               {"x": 5.5, "y": 9.6, "phi": 90}])
            again = Project.open(project.root)
            self.assertAlmostEqual(again.unit_cell["cell_x"], 11.1)
            self.assertEqual(len(again.unit_cell["molecules"]), 2)
            again.clear_unit_cell()
            self.assertIsNone(again.unit_cell)

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
